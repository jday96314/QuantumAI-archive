"""Use remaining tensor-error budget to lower the optimized SOS gap.

The gap's envelope gradient is a rank-one tensor in the projected occupied
one-body subspace. A positive update preserves the PSD physical factors.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--target',type=float,default=3.916);p.add_argument('--cap',type=float,default=.01999)
p.add_argument('--cutoff',type=float,default=1e-8);a=p.parse_args();HERE=Path(__file__).resolve().parent
z=dict(np.load(a.fit));ref=np.load(HERE/'parameters/chem_reference.npz');sp=np.load(HERE/'parameters/reference_spectrum.npz')
u,w=z['U'][0],z['w'][0];B,C=w.shape;assert z['w'].shape[0]==1
ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];P=u[:,ii]*u[:,jj]*fac
X=w.T@P;S=X@X.T;G=(sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T
res=X.T@X-G;err0=np.linalg.norm(res)
xf=np.einsum('bp,bq,bc->cpq',u,u,w,optimize=True).reshape(C,-1)
h0=ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(76)
constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
off=z['W'][0];h=h0-(off@xf).reshape(76,76);e,v=eigh(h)
Z=(v*(np.sign(e)-1))@v.T;vproj=X@(Z[ii,jj]*fac)
ee,vv=eigh(S);q=(vv*(1/np.maximum(ee,a.cutoff)))@(vv.T@vproj)
Y=q@X;yn=float(Y@Y);q2=float(q@q);cross=float(Y@res@Y)
maxeps=(-cross+np.sqrt(max(0,cross*cross+yn*yn*(a.cap*a.cap-err0*err0))))/(yn*yn)
print(json.dumps(dict(error_2=float(err0),direction_norm=yn,coefficient_norm=q2,max_epsilon=maxeps,
                      max_linear_gap_reduction=.5*maxeps*yn*yn)),flush=True)
def trial(epsilon):
 m=np.eye(C)+(np.sqrt(1+epsilon*q2)-1)/q2*np.outer(q,q)
 nw=w@m;nxf=m.T@xf;tr=np.trace(nxf.reshape(C,76,76),axis1=1,axis2=2)
 def fun(o,grad=False):
  hh=h0-(o@nxf).reshape(76,76);ev,ve=eigh(.5*(hh+hh.T));r=np.sqrt(ev*ev+1e-14)
  gap=constant+r.sum()-ev.sum()+.5*o@o
  if grad:return gap,o+tr-nxf@((ve*(ev/r))@ve.T).reshape(-1)
  lam=r.sum()+.25*np.square(abs(nw).sum(0)+abs(o+nw.sum(0))).sum()
  return gap,lam
 init=np.linalg.solve(m,off)
 rr=minimize(lambda x:fun(x,True),init,jac=True,method='L-BFGS-B',options=dict(maxiter=2000,ftol=1e-14,gtol=1e-8))
 gap,lam=fun(rr.x)
 err=np.sqrt(max(0,err0**2+2*epsilon*cross+epsilon*epsilon*yn*yn))
 return gap,lam,err,nw,rr.x
lo,hi=0.,maxeps;best=None
for step in range(14):
 eps=hi if step==0 else .5*(lo+hi);g,l,err,nw,offnew=trial(eps)
 print(json.dumps(dict(step=step,epsilon=eps,error_2=float(err),gap=float(g),Lambda=float(l),
                       lambda_eff=float(np.sqrt(max(0,g*(2*l-g)))))),flush=True)
 if g>=a.target:lo=eps;best=(nw,offnew,err,g)
 else:hi=eps
 if step==0 and g>a.target:break
if best is None:
 g,l,err,nw,offnew=trial(lo);best=(nw,offnew,err,g)
z['w']=best[0][None,:,:];z['W']=best[1][None,:];z['err2']=best[2]
np.savez(a.out,**z)
