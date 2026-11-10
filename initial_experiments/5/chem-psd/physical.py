"""Exact decoded refits, offset/gap optimization and native word evaluation.

Offset formula follows completed chem-qfit/physical.py and the unchanged host.
Fixed S and beta_1 are always taken from load_inputs, never optimized.
"""
import argparse,time,subprocess
from dataclasses import replace
from scipy.linalg import eigh
from scipy.optimize import minimize,brentq
from core import *
from solve import solve
from angle_gauge_probe import encode
from dfthcsos.verify import evaluate_words,passes_accuracy,resource_costs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x,build_h1_prime
from signed_oracle import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from frontier import price

def offset_objective(off,xf,tr,h0,constant,w,gradient=True):
    h=h0-(off@xf).reshape(76,76);ev,V=eigh((h+h.T)/2)
    ab=np.sqrt(ev*ev+1e-16);gap=constant+ab.sum()-ev.sum()+.5*(off@off)
    if gradient:return gap,off+tr-xf@((V*(ev/ab))@V.T).ravel()
    lam=ab.sum()+.25*((np.abs(w).sum(1).ravel()+np.abs(off+w.sum(1).ravel()))**2).sum()
    return float(gap),float(lam)

@torch.no_grad()
def onebody_refine(h,tree,steps=1000,target=.015):
    """Exact integer +/-1 proposals on the weighted one-body matrix.

    Algebra adapted from completed chem-qfit and round-3 quantize_onebody;
    proposals are accepted using full matrix reconstruction, not angle loss.
    """
    ev,V=eigh((h+h.T)/2);theta=tensor(encode(V.T,tree,8,True));e=tensor(ev);H=tensor(h);history=[]
    alpha=np.pi/128
    for step in range(steps):
        v=torch.zeros((76,76),**D);v[:,0]=1;J=torch.zeros((76,75,76),**D);Q=J.clone()
        for t,(p,c) in enumerate(tree):
            co=(theta[:,t]*alpha).cos();sn=(theta[:,t]*alpha).sin()
            for derivative in [J,Q]:
                pp,cc=derivative[:,:,p].clone(),derivative[:,:,c].clone()
                derivative[:,:,p]=co[:,None]*pp-sn[:,None]*cc;derivative[:,:,c]=sn[:,None]*pp+co[:,None]*cc
            pp,cc=v[:,p].clone(),v[:,c].clone();v[:,p]=co*pp-sn*cc;v[:,c]=sn*pp+co*cc
            J[:,t,p]=-v[:,c];J[:,t,c]=v[:,p];Q[:,t,p]=-v[:,p];Q[:,t,c]=-v[:,c]
        residual=(v.T*e)@v-H;loss=float(residual.square().sum());history.append(loss**.5)
        if loss<=target*target:break
        can=v[:,None,None,:]+np.sin(alpha)*torch.stack([J,-J],2)+(1-np.cos(alpha))*Q[:,:,None,:]
        old=(v*(v@residual)).sum(1);quad=(can*(can@residual)).sum(3);dot=(can*v[:,None,None,:]).sum(3)
        delta=2*e[:,None,None]*(quad-old[:,None,None])+2*e[:,None,None].square()*(1-dot.square())
        value,index=delta.reshape(76,-1).min(1);order=value.argsort();order=order[value[order]<-1e-14];changed=False
        for count in sorted(set([len(order),32,16,8,4,2,1]),reverse=True):
            if not 0<count<=len(order):continue
            rows=order[:count];cols=index[rows]//2;signs=index[rows]%2;trial=v.clone();trial[rows]=can[rows,cols,signs]
            if float(((trial.T*e)@trial-H).square().sum())<loss-1e-14:
                theta[rows,cols]+=torch.where(signs==0,1.,-1.);changed=True;break
        if not changed:break
    return encode(reconstruct_u_vectors(theta.cpu().numpy()/128,tree),tree,8,True),history

def optimize_offsets(inp,u,w,tag,jacobi=0,task=False):
    X=reconstruct_x(u,w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).ravel()
    h0=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],np.zeros(w.shape[::2]),X)
    const=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1
    fun=lambda o:offset_objective(o,xf,tr,h0,const,w)
    res=minimize(fun,np.zeros(xf.shape[0]),jac=True,method='L-BFGS-B',options=dict(maxiter=1000,ftol=1e-14,gtol=1e-9,maxls=40))
    off=res.x;mingap=fun(off)[0]
    # Enforce a small strict gap-floor margin along an explicit feasible line.
    if mingap<3.9152:
        origin=off.copy();direction=-w.sum(1).ravel()-origin
        if fun(origin+direction)[0]<3.9152:direction=-origin
        length=1.
        while fun(origin+length*direction)[0]<3.9152:length*=2
        t=brentq(lambda t:fun(origin+t*direction)[0]-3.9152,0,length,xtol=1e-13)
        off=origin+t*direction
    off=off.reshape(w.shape[0],w.shape[2]);checks=[]
    if jacobi:
        from jacobi_l1 import encode as write_binary,decode as read_binary
        for r in range(w.shape[0]):
            A=np.r_[w[r],(off[r]+w[r].sum(0))[None]]
            src=HERE/f'parameters/{tag}_gauge{r}_in.bin';dst=HERE/f'parameters/{tag}_gauge{r}_out.bin'
            write_binary(A,src)
            with (HERE/f'results/{tag}_gauge{r}.log').open('w') as log:
                subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),str(jacobi),'12000','52719'],stdout=log,check=True)
            rot=read_binary(dst);err=float(abs(rot@rot.T-A@A.T).max());assert err<1e-8
            checks.append(dict(gram_max_error=err,lambda_h2_before=float(.25*(abs(A).sum(0)**2).sum()),lambda_h2_after=float(.25*(abs(rot).sum(0)**2).sum())))
            w[r]=rot[:-1];off[r]=rot[-1]-rot[:-1].sum(0)
    task_result=None
    if task:
        X=reconstruct_x(u,w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).ravel()
        def task_values(o):
            h=h0-(o@xf).reshape(76,76);ev,V=eigh((h+h.T)/2);ab=np.sqrt(ev*ev+1e-16)
            dnorm=-xf@((V*(ev/ab))@V.T).ravel();b=o+w.sum(1).ravel();bw=np.sqrt(b*b+1e-16)
            gap=const+ab.sum()-ev.sum()+.5*(o@o);lam=ab.sum()+.25*((abs(w).sum(1).ravel()+bw)**2).sum()
            dg=o+tr+dnorm;dl=dnorm+.5*(abs(w).sum(1).ravel()+bw)*b/bw
            return gap,lam,dg,dl
        def objective(o):
            g,l,dg,dl=task_values(o);return g*(2*l-g),2*(l-g)*dg+2*g*dl
        cons=[dict(type='ineq',fun=lambda o:task_values(o)[0]-3.9152,jac=lambda o:task_values(o)[2]),
              dict(type='ineq',fun=lambda o:task_values(o)[1]-49.00001,jac=lambda o:task_values(o)[3]),
              dict(type='ineq',fun=lambda o:2*task_values(o)[1]-task_values(o)[0]-1e-7,jac=lambda o:2*task_values(o)[3]-task_values(o)[2])]
        before=objective(off.ravel())[0]
        opt=minimize(objective,off.ravel(),jac=True,method='SLSQP',constraints=cons,options=dict(maxiter=400,ftol=1e-9))
        g,l,_,_=task_values(opt.x)
        if opt.success and g>=3.9152-1e-9 and l>49 and 0<g<2*l and objective(opt.x)[0]<before:off=opt.x.reshape(off.shape)
        task_result=dict(termination=str(opt.message),radicand_before=float(before),radicand_after=float(objective(off.ravel())[0]))
    return w,off,dict(minimum_gap=mingap,offset_termination=str(res.message),gauge=checks,task_optimization=task_result)

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--tag',required=True);p.add_argument('--rank',type=int,required=True)
    p.add_argument('--steps',type=int,default=1200);p.add_argument('--jacobi',type=int,default=0);p.add_argument('--build',action='store_true');a=p.parse_args()
    z=dict(np.load(a.input));u=z['U'];R,B,_=u.shape;tree=z['tree_edges'];C=a.rank
    m=Covariance(u);S=m.from_weights(z['w']) if 'w' in z else None
    S,hc=solve(m,S,a.steps,C,'coordinate',200)
    wc=m.weights(S,C).cpu().numpy();continuous=m.diagnostics(S,True,C)
    np.savez(HERE/f'parameters/{a.tag}_continuous.npz',U=u,w=wc,S=S.cpu().numpy(),tree_edges=tree)
    words=encode(u.reshape(-1,76),tree,8,True);uq=reconstruct_u_vectors(words/128,tree).reshape(R,B,76)
    mq=Covariance(uq);Sq=mq.from_weights(wc);rounded=mq.diagnostics(Sq,True)
    Sq,hq=solve(mq,Sq,a.steps,C,'coordinate',200)
    w=mq.weights(Sq,C).cpu().numpy();quantized=mq.diagnostics(Sq,True,C)
    np.savez(HERE/f'parameters/{a.tag}_quantized.npz',U=uq,w=w,S=Sq.cpu().numpy(),tree_edges=tree,angle_words=words)
    inp=load_inputs(INC);base=dict(np.load(INC));w,off,offsets=optimize_offsets(inp,uq,w,a.tag,a.jacobi,quantized['direct_error_2']<=.02)
    X=reconstruct_x(uq,w);h=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],off,X)
    _,V=eigh((h+h.T)/2);one=encode(V.T,tree,8,True)
    allwords=np.r_[one,words]
    params=dict(U=uq,W=w,H2_I=off+w.sum(1),H1_I=base['H1_I'],Hsym_lower=base['Hsym_lower'])
    this=replace(inp,params=params,w=w,W=off);metrics=evaluate_words(this,allwords,tree,8)
    onebody_history=[]
    if metrics['error_2']<=.02 and metrics['error_1']>.015:
        one,onebody_history=onebody_refine(h,tree);allwords=np.r_[one,words]
        metrics=evaluate_words(this,allwords,tree,8)
    assert abs(metrics['error_2']-quantized['direct_error_2'])<1e-7
    np.savez(HERE/f'parameters/{a.tag}_physical.npz',**params,angle_words=allwords,tree_edges=tree)
    resources=price(R,B,C,8,7);cost=metrics['lambda_eff']*resources['T_total_bound']*resources['Q_total_bound']
    rec=dict(input=a.input,rank=C,continuous=continuous,rounded=rounded,quantized=quantized,offsets=offsets,
       native_metrics=metrics,onebody_history=onebody_history,passes_physical_gates=passes_accuracy(metrics),resource_bound=resources,
       conditional_invalid_product=cost if not passes_accuracy(metrics) else None,
       error_margins=dict(error_1=.02-metrics['error_1'],error_2=.02-metrics['error_2']),accepted_task_score=None)
    if a.build and passes_accuracy(metrics):
        out=HERE/f'candidates/{a.tag}';out.mkdir(parents=True,exist_ok=True)
        np.savez(out/'submission.npz',**params,angle_words=allwords,tree_edges=tree)
        circuits=[]
        for k,delta in enumerate(batch_deltas(allwords,7)):
            circuit=build(delta,8,resources['lowbits']);assert circuit.num_qubits<=128
            circuit.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');circuits.append(circuit)
        T=2*sum(c.operation_statistics().num_toffoli_like for c in circuits)
        rec['compiled_resources']=resource_costs(this,circuits,8,T)
        rec['native_verification_pending']=str(out)
    report(HERE/f'results/{a.tag}_physical.json',rec)
if __name__=='__main__':main()
