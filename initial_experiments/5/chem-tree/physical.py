"""Native normalization, fixed-shift physics and actual complete resource screens.

Offset and one-body algorithms adapt completed chem-qfit/physical.py. Native
load_inputs/evaluate_words/resource_costs remain unchanged. The read-only
round-4 Jacobi binary is used only for physical-Gram-preserving column gauges.
"""
import argparse,subprocess
from dataclasses import replace
from scipy.linalg import eigh
from scipy.optimize import minimize
from core import *
from phase_grid import phase_encode
from dfthcsos.verify import evaluate_words,passes_accuracy,resource_costs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x,build_h1_prime
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
from oracles import build

def offset_fun(off,xf,tr,h0,constant,w,grad=True):
    h=h0-(off@xf).reshape(N,N);ev,V=eigh((h+h.T)/2)
    sm=np.sqrt(ev*ev+1e-14);gap=constant+sm.sum()-ev.sum()+.5*np.dot(off,off)
    if grad:return gap,off+tr-xf@((V*(ev/sm))@V.T).reshape(-1)
    lam=sm.sum()+.25*((np.abs(w).sum(1).reshape(-1)+np.abs(off+w.sum(1).reshape(-1)))**2).sum()
    return gap,lam

def onebody(h,tree,bits=8,phase=None,steps=300):
    ev,V=np.linalg.eigh((h+h.T)/2);stride=1<<(bits-8)
    if phase is None:phase=np.zeros(75,dtype=np.int64)
    words,u=phase_encode(V.T,tree,phase,bits)
    initial=None;start=time.time();dt=np.pi/128
    for step in range(steps):
        v=np.zeros((N,N));v[:,0]=1;J=np.zeros((N,75,N));Q=J.copy()
        for j,(p,c) in enumerate(tree):
            co=np.cos(words[:,j]*np.pi/(1<<(bits-1)));sn=np.sin(words[:,j]*np.pi/(1<<(bits-1)))
            for Dm in [J,Q]:
                a=Dm[:,:,p].copy();b=Dm[:,:,c].copy();Dm[:,:,p]=co[:,None]*a-sn[:,None]*b;Dm[:,:,c]=sn[:,None]*a+co[:,None]*b
            a=v[:,p].copy();b=v[:,c].copy();v[:,p]=co*a-sn*b;v[:,c]=sn*a+co*b
            J[:,j,p]=-v[:,c];J[:,j,c]=v[:,p];Q[:,j,p]=-v[:,p];Q[:,j,c]=-v[:,c]
        E=(v.T*ev)@v-h;err=float(np.linalg.norm(E))
        if initial is None:initial=err
        if err<.015:break
        can=v[:,None,None,:]+np.sin(dt)*np.stack([J,-J],2)+(1-np.cos(dt))*Q[:,:,None,:]
        old=np.sum(v*(v@E),axis=-1);quad=np.sum(can*(can@E),axis=-1);dot=np.sum(can*v[:,None,None,:],axis=-1)
        delta=2*ev[:,None,None]*(quad-old[:,None,None])+2*ev[:,None,None]**2*(1-dot**2)
        index=delta.reshape(N,-1).argmin(1);value=delta.reshape(N,-1).min(1)
        order=value.argsort();order=order[value[order]<-1e-14];changed=False
        for count in sorted(set([len(order),32,16,8,4,2,1]),reverse=True):
            if count<=0 or count>len(order):continue
            rr=order[:count];cc=index[rr]//2;ss=index[rr]%2;new=v.copy();new[rr]=can[rr,cc,ss]
            if np.linalg.norm((new.T*ev)@new-h)**2<err**2-1e-14:
                words[rr,cc]+=(np.where(ss==0,1,-1)*stride);words%=1<<bits;changed=True;break
        if not changed:break
    words,_=phase_encode(decode(words*np.pi/(1<<(bits-1)),tree),tree,phase,bits)
    return words,dict(initial=initial,final=err,steps=step,seconds=time.time()-start)

def assess(fit,jacobi=8,build_circuits=True):
    z,u,w,tree,words=load(fit);R,B,C=w.shape;bits=int(z.get('brot',8));phase=z.get('grid_phase',np.zeros(75,dtype=np.int64))
    if words is None:words=phase_encode(u,tree,phase,bits)[0].reshape(-1,75)
    u=decode(words.reshape(R,B,75)*np.pi/(1<<(bits-1)),tree)
    old=u.copy();words,u=phase_encode(u,tree,phase,bits);words=words.reshape(-1,75)
    assert np.max(np.abs(np.abs(np.sum(old*u,axis=-1))-1))<1e-12
    inp=native();base=np.load(INC);start=time.time()
    h0=inp.ref['h1_thc']-.5*inp.ref['electrons']*inp.S+inp.beta_1*np.eye(N)
    constant=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1
    X=reconstruct_x(u,w);xf=X.reshape(-1,N*N);tr=np.trace(X,axis1=-2,axis2=-1).reshape(-1)
    fun=lambda off:offset_fun(off,xf,tr,h0,constant,w)
    res=minimize(fun,np.zeros(R*C),jac=True,method='L-BFGS-B',options=dict(maxiter=800,ftol=1e-14,gtol=1e-9,maxls=40))
    off=res.x;minimum_gap=float(fun(off)[0])
    if minimum_gap<3.9152:
        origin=off.copy();direction=-w.sum(1).reshape(-1)-origin
        if fun(origin+direction)[0]<3.9152:direction=-origin
        lo,hi=0.,1.
        for _ in range(70):
            mid=(lo+hi)/2
            if fun(origin+mid*direction)[0]<3.9152:lo=mid
            else:hi=mid
        off=origin+hi*direction
    off=off.reshape(R,C);h=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],off,X)
    gauge=[];stem=Path(fit).stem
    if jacobi:
        from jacobi_l1 import encode as bin_encode,decode as bin_decode
        for r in range(R):
            A=np.concatenate([w[r],(off[r]+w[r].sum(0))[None]],axis=0)
            src=HERE/f'candidates/{stem}_gauge{r}_input.bin';dst=src.with_name(src.name.replace('input','output'));bin_encode(A,src)
            with (HERE/f'results/{stem}_gauge{r}.log').open('w') as log:
                subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),str(jacobi),'12000','791291'],stdout=log,check=True)
            rot=bin_decode(dst);delta=rot@rot.T-A@A.T;assert abs(delta).max()<1e-9
            gauge.append(dict(gram_max_error=float(abs(delta).max()),lambda_h2_before=float(.25*(np.abs(A).sum(0)**2).sum()),lambda_h2_after=float(.25*(np.abs(rot).sum(0)**2).sum())))
            w[r]=rot[:-1];off[r]=rot[-1]-rot[:-1].sum(0)
    one,oc=onebody(h,tree,bits,phase)
    allwords=np.concatenate([one,words.reshape(-1,75)])
    params=dict(U=u,W=w,H2_I=off+w.sum(1),H1_I=base['H1_I'],Hsym_lower=base['Hsym_lower'])
    this=replace(inp,params=params,w=w,W=off)
    metrics=evaluate_words(this,allwords,tree,bits)
    record=dict(fit=str(fit),shape=[R,B,C],bits=bits,metrics=metrics,physical_pass=passes_accuracy(metrics),
                minimum_gap=minimum_gap,offset_success=bool(res.success),offset_message=res.message,onebody=oc,gauge=gauge)
    if build_circuits:
        out=HERE/'candidates'/f'{stem}_submission';out.mkdir(exist_ok=True)
        np.savez(out/'submission.npz',**params,tree_edges=tree,angle_words=allwords,brot=bits,grid_phase=phase)
        circuits=[];J=7 if bits==8 else 6
        records=[]
        for k,delta in enumerate(batch_deltas(allwords,J)):
            c=build(delta,bits,bits-8);assert c.num_qubits<=128
            c.to_file(str(out/f'angle_k_{k}.kmb'),format='kmb');circuits.append(c)
            records.append(dict(k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
        T=2*sum(r['T'] for r in records)
        record['resources']=resource_costs(this,circuits,bits,T)
        record['conditional_product']=metrics['lambda_eff']*record['resources']['toffoli']*record['resources']['qubits']
        record['lambda_break_even']=708547388.9553173/record['resources']['toffoli']/record['resources']['qubits']
        record['oracles']=records;record['submission_directory']=str(out)
    record['seconds']=time.time()-start
    # This is a physical/resource screen, never a native PASS assertion.
    report(HERE/f'results/{stem}_physical.json',record)
    np.savez(HERE/f'candidates/{stem}_physical.npz',U=u,w=w,W=off,tree_edges=tree,angle_words=words.reshape(-1,75),onebody_words=one,brot=bits,grid_phase=phase)
    return record

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('fits',nargs='+');p.add_argument('--jacobi',type=int,default=8);p.add_argument('--no-build',action='store_true');a=p.parse_args()
    for fit in a.fits:assess(fit,a.jacobi,not a.no_build)
