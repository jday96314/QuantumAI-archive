"""Analytic CPU optimization for the small fixed-factor SOS gap problem."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import eigh
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
a=p.parse_args();HERE=Path(__file__).resolve().parent;z=dict(np.load(a.fit));ref=np.load(HERE/'parameters/chem_reference.npz')
u,w=z['U'],z['w'];R,B,C=w.shape;n=76
X=np.einsum('rbp,rbq,rbc->rcpq',u,u,w,optimize=True).reshape(-1,n,n);xf=X.reshape(R*C,-1)
h0=ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(n)
constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
traces=np.trace(X,axis1=1,axis2=2)
def vals(off,gradient=False):
    h=h0-(off@xf).reshape(n,n);e,v=eigh(.5*(h+h.T));r=np.sqrt(e*e+1e-14)
    gap=constant+r.sum()-e.sum()+.5*np.square(off).sum()
    if gradient:
        sign=(v*(e/r))@v.T
        return gap,off+traces-xf@sign.reshape(-1)
    wb=off.reshape(R,C)+w.sum(1);lam=r.sum()+.25*np.square(abs(w).sum(1)+abs(wb)).sum()
    return gap,lam
start=time.time();res=minimize(lambda x:vals(x,True),np.zeros(R*C),jac=True,method='L-BFGS-B',
                             options=dict(maxiter=2000,ftol=1e-14,gtol=1e-8,maxls=40))
W=res.x;g,l=vals(W)
if g<3.916:
    origin=W.copy();direction=-w.sum(1).reshape(-1)-origin
    if vals(origin+direction)[0]<3.916:direction=-origin
    lo,hi=0.,1.
    for _ in range(80):
        mid=(lo+hi)/2
        if vals(origin+mid*direction)[0]<3.916:lo=mid
        else:hi=mid
    W=origin+hi*direction
g,l=vals(W);z['W']=W.reshape(R,C);np.savez(a.out,**z)
A=np.concatenate([w,(z['W']+w.sum(1))[:,None,:]],1);gram=A@A.transpose(0,2,1)
print(json.dumps(dict(gap=float(g),Lambda=float(l),lambda_eff=float(np.sqrt(g*(2*l-g))),
                     gram_lower_bound=float(.25*abs(gram).sum()),iterations=res.nit,success=bool(res.success),seconds=time.time()-start)),flush=True)
