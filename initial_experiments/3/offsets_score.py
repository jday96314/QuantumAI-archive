"""Optimize the actual amplified normalization with fixed physical factors."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigh
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--starts',type=int,default=4);a=p.parse_args();HERE=Path(__file__).resolve().parent
z=dict(np.load(a.fit));ref=np.load(HERE/'parameters/chem_reference.npz');u,w=z['U'],z['w'];R,B,C=w.shape
X=np.einsum('rbp,rbq,rbc->rcpq',u,u,w,optimize=True).reshape(-1,76,76);xf=X.reshape(R*C,-1)
h0=ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(76)
constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
traces=np.trace(X,axis1=1,axis2=2);s=w.sum(1).reshape(-1);l1=abs(w).sum(1).reshape(-1)
def vals(off):
 h=h0-(off@xf).reshape(76,76);e,v=eigh(.5*(h+h.T));r=np.sqrt(e*e+1e-14)
 sign=(v*(e/r))@v.T;xs=xf@sign.reshape(-1)
 gap=constant+r.sum()-e.sum()+.5*np.square(off).sum();dg=off+traces-xs
 c=off+s;lam=r.sum()+.25*np.square(l1+abs(c)).sum();dl=-xs+.5*np.sign(c)*(l1+abs(c))
 score=gap*(2*lam-gap);ds=2*gap*dl+2*(lam-gap)*dg
 return score,ds,gap,dg,lam
start=time.time();best=vals(z['W'].reshape(-1))[0];bestw=z['W'].reshape(-1).copy()
cons={'type':'ineq','fun':lambda x:vals(x)[2]-3.9155,'jac':lambda x:vals(x)[3]}
rng=np.random.default_rng(99)
for k in range(a.starts):
 init=bestw.copy() if k==0 else (-s if k==1 else bestw+rng.normal(size=len(s))*.25)
 res=minimize(lambda x:vals(x)[:2],init,jac=True,method='SLSQP',constraints=cons,
              options=dict(maxiter=600,ftol=1e-8))
 val,_,gap,_,lam=vals(res.x)
 print(json.dumps(dict(start=k,success=bool(res.success),steps=res.nit,gap=float(gap),Lambda=float(lam),
                       lambda_eff=float(np.sqrt(max(0,val))),seconds=time.time()-start)),flush=True)
 if val<best and gap>=3.91549 and val>0:best=val;bestw=res.x.copy()
z['W']=bestw.reshape(R,C);np.savez(a.out,**z)
