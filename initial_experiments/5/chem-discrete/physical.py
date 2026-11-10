"""Native physical assessment, exact-semantic offset solver and ANGLE builder.

Offset objective/gradient adapted from AGY fast_optimize_offsets.py (round 4),
generalized to all R blocks and using native fixed shifts and reconstructed X.
One-body finite rotation updates follow experiment 3 quantize_onebody.py.
"""
import argparse,json,time,subprocess
from scipy.linalg import eigh
from scipy.optimize import minimize
from core import *
from oracles import build
from dfthcsos.verify import load_inputs,evaluate_words,passes_accuracy,resource_costs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x,build_h1_prime
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas

def offset_fun(off,xf,tr,h0,constant,w,grad=True):
    h=h0-(off@xf).reshape(76,76);ev,V=eigh((h+h.T)/2)
    sm=np.sqrt(ev*ev+1e-14);gap=constant+sm.sum()-ev.sum()+.5*np.dot(off,off)
    if grad:return gap,off+tr-xf@((V*(ev/sm))@V.T).reshape(-1)
    lam=sm.sum()+.25*((np.abs(w).sum(1).reshape(-1)+np.abs(off+w.sum(1).reshape(-1)))**2).sum()
    return gap,lam

def onebody(h,tree,steps=800,target=.015):
    from types import SimpleNamespace
    ev,V=np.linalg.eigh((h+h.T)/2);words=encode(V.T,tree,8,True)
    # These 76-row operations are launch-bound on the GPU and run faster on
    # the CPU while large dictionary searches occupy the device.
    d=dict(device='cpu',dtype=torch.float64);m=SimpleNamespace(d=d,tensor=lambda x:torch.tensor(x,**d))
    t=m.tensor(words);e=m.tensor(ev);H=m.tensor(h);start=time.time();initial=None
    with torch.no_grad():
        for step in range(steps):
            v=torch.zeros((76,76),**m.d);v[:,0]=1;J=torch.zeros((76,75,76),**m.d);Q=J.clone()
            for j,(p,c) in enumerate(tree):
                co=(t[:,j]*np.pi/128).cos();sn=(t[:,j]*np.pi/128).sin()
                for D in [J,Q]:
                    pp,cc=D[:,:,p].clone(),D[:,:,c].clone();D[:,:,p]=co[:,None]*pp-sn[:,None]*cc;D[:,:,c]=sn[:,None]*pp+co[:,None]*cc
                pp,cc=v[:,p].clone(),v[:,c].clone();v[:,p]=co*pp-sn*cc;v[:,c]=sn*pp+co*cc
                J[:,j,p]=-v[:,c];J[:,j,c]=v[:,p];Q[:,j,p]=-v[:,p];Q[:,j,c]=-v[:,c]
            D=(v.T*e)@v-H;err=float(D.norm())
            if initial is None:initial=err
            if err<=target:break
            can=v[:,None,None,:]+np.sin(np.pi/128)*torch.stack([J,-J],2)+(1-np.cos(np.pi/128))*Q[:,:,None,:]
            old=(v*(v@D)).sum(-1);quad=(can*(can@D)).sum(-1);dot=(can*v[:,None,None,:]).sum(-1)
            delta=2*e[:,None,None]*(quad-old[:,None,None])+2*e[:,None,None].square()*(1-dot.square())
            value,ix=delta.reshape(76,-1).min(1);order=value.argsort();order=order[value[order]<-1e-14];changed=False
            for count in sorted(set([len(order),32,16,8,4,2,1]),reverse=True):
                if count<=0 or count>len(order):continue
                rr=order[:count];cc=ix[rr]//2;ss=ix[rr]%2;new=v.clone();new[rr]=can[rr,cc,ss]
                if float(((new.T*e)@new-H).square().sum())<err*err-1e-14:
                    t[rr,cc]+=torch.where(ss==0,1.,-1.);changed=True;break
            if not changed:break
    # Re-gauge to the oracle's signed extension without projector changes.
    out=encode(reconstruct_u_vectors(t.cpu().numpy()/128,tree),tree,8,True)
    return out,dict(initial=initial,final=err,steps=step,seconds=time.time()-start)

def main():
    p=argparse.ArgumentParser();p.add_argument('fits',nargs='+');p.add_argument('--build',action='store_true');p.add_argument('--resources',action='store_true');p.add_argument('--keep-onebody',action='store_true');p.add_argument('--skip-offset',action='store_true');p.add_argument('--preserve-offset',action='store_true');p.add_argument('--align',action='store_true');p.add_argument('--jacobi',type=int,default=0);a=p.parse_args()
    inp=load_inputs(INC/'submission.npz');base=dict(np.load(INC/'submission.npz'))
    h0=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],np.zeros_like(inp.W),reconstruct_x(base['U'],inp.w))
    constant=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1
    # Independent check of AGY's exact objective, gradient and native semantics.
    X=reconstruct_x(base['U'],inp.w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).reshape(-1);off=inp.W.reshape(-1)
    gap,lam=offset_fun(off,xf,tr,h0,constant,inp.w,False);_,grad=offset_fun(off,xf,tr,h0,constant,inp.w)
    rng=np.random.default_rng(541);checks=[]
    for _ in range(4):
        v=rng.normal(size=len(off));v/=np.linalg.norm(v);eps=1e-5
        fd=(offset_fun(off+eps*v,xf,tr,h0,constant,inp.w)[0]-offset_fun(off-eps*v,xf,tr,h0,constant,inp.w)[0])/(2*eps)
        checks.append(dict(finite_difference=fd,analytic=float(grad@v),abs_error=abs(fd-grad@v)))
    native=evaluate_words(inp,base['angle_words'],base['tree_edges'],8)
    check=dict(gap=gap,Lambda=lam,native_gap=native['E_gap'],native_Lambda=native['Lambda'],gradient_checks=checks)
    assert abs(gap-native['E_gap'])<1e-8 and abs(lam-native['Lambda'])<1e-8
    (HERE/'results/offset_semantics.json').write_text(json.dumps(check,indent=2)+'\n')
    for fit in a.fits:
        z,u,w,tree,words=load(fit);R,B,C=w.shape
        if words is None:words=encode(u.reshape(-1,76),tree,8,True)
        u=reconstruct_u_vectors(words/128,tree).reshape(R,B,76)
        words=encode(u.reshape(-1,76),tree,8,True);u=reconstruct_u_vectors(words/128,tree).reshape(R,B,76)
        if a.align and R==1:
            nearest=np.abs(u[0]@base['U'][0].T).argmax(1)
            warm=np.zeros((B,C));warm[:,:min(C,inp.w.shape[-1])]=inp.w[0,nearest,:min(C,inp.w.shape[-1])]
            left,_,right=np.linalg.svd(w[0].T@warm);w=w@(left@right)
        X=reconstruct_x(u,w);xf=X.reshape(-1,76*76);tr=np.trace(X,axis1=-2,axis2=-1).reshape(-1)
        if a.preserve_offset:off=np.asarray(z['W']).reshape(-1).copy()
        elif a.skip_offset and R==1 and C==350:off=inp.W.reshape(-1).copy()
        else:
            fun=lambda o:offset_fun(o,xf,tr,h0,constant,w)
            res=minimize(fun,np.zeros(R*C),jac=True,method='L-BFGS-B',options=dict(maxiter=800,ftol=1e-14,gtol=1e-9,maxls=40))
            off=res.x
            # Move on a concrete line to a small strict-floor margin. This is
            # needed if the convex gap minimizer lies below the accepted floor.
            if fun(off)[0]<3.9152:
                origin=off.copy();direction=-w.sum(1).reshape(-1)-origin
                if fun(origin+direction)[0]<3.9152:direction=-origin
                lo,hi=0.,1.
                for _ in range(70):
                    mid=(lo+hi)/2
                    if fun(origin+mid*direction)[0]<3.9152:lo=mid
                    else:hi=mid
                off=origin+hi*direction
        off=off.reshape(R,C);h=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],off,X)
        gauge=[]
        if a.jacobi:
            from jacobi_l1 import encode as bin_encode,decode as bin_decode
            for r in range(R):
                A=np.concatenate([w[r],(off[r]+w[r].sum(0))[None]],axis=0)
                src=HERE/f'candidates/{Path(fit).stem}_gauge{r}_input.bin';dst=src.with_name(src.name.replace('input','output'))
                bin_encode(A,src)
                with (HERE/f'results/{Path(fit).stem}_gauge{r}.log').open('w') as log:
                    subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),str(a.jacobi),'12000','791291'],stdout=log,check=True)
                rot=bin_decode(dst);delta=rot@rot.T-A@A.T;assert abs(delta).max()<1e-9
                gauge.append(dict(gram_max_error=float(abs(delta).max()),lambda_h2_before=float(.25*(np.abs(A).sum(0)**2).sum()),lambda_h2_after=float(.25*(np.abs(rot).sum(0)**2).sum())))
                w[r]=rot[:-1];off[r]=rot[-1]-rot[:-1].sum(0)
        if a.keep_onebody and 'onebody_angle_words' in z:
            one=z['onebody_angle_words'];oc=dict(preserved=True)
        else:one,oc=onebody(h,tree)
        allwords=np.concatenate([one,words]);params=dict(U=u,W=w,H2_I=off+w.sum(1),H1_I=base['H1_I'],Hsym_lower=base['Hsym_lower'])
        # Reuse the already loaded, fixed-shift native reference.
        from dataclasses import replace
        this=replace(inp,params=params,w=w,W=off)
        metrics=evaluate_words(this,allwords,tree,8)
        record=dict(fit=fit,shape=[R,B,C],metrics=metrics,physical_pass=passes_accuracy(metrics),onebody=oc,gauge=gauge)
        assessed=HERE/'candidates'/f'evaluated_{Path(fit).stem}.npz'
        save(assessed,u,w,tree,words,W=off,onebody_angle_words=one,err2=metrics['error_2'])
        record['evaluated_checkpoint']=str(assessed)
        record['native_status']='not_tested' if record['physical_pass'] else 'physical_reject'
        if a.resources or (a.build and record['physical_pass']):
            out=HERE/'candidates'/Path(fit).stem
            if a.build and record['physical_pass']:
                out.mkdir(exist_ok=True)
                np.savez(out/'submission.npz',**params,tree_edges=tree,angle_words=allwords)
            circuits=[]
            for k,delta in enumerate(batch_deltas(allwords,7)):
                circuit=build(delta);assert circuit.num_qubits<=128
                if a.build and record['physical_pass']:circuit.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb')
                circuits.append(circuit)
            T=2*sum(c.operation_statistics().num_toffoli_like for c in circuits)
            record['resources_static']=resource_costs(this,circuits,8,T)
            record['estimated_score']=metrics['lambda_eff']*record['resources_static']['toffoli']*record['resources_static']['qubits']
            record['break_even_lambda']=708547388.9553/record['resources_static']['toffoli']/record['resources_static']['qubits']
            if a.build and record['physical_pass']:record['circuit_directory']=str(out)
        print(json.dumps(record),flush=True)
        (HERE/f'results/{Path(fit).stem}_physical.json').write_text(json.dumps(record,indent=2)+'\n')

if __name__=='__main__':main()
