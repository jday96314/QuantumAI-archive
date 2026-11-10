"""Refit nonnegative ensemble strengths with decoded angle vectors fixed.

The two-body tensor is linear in squared column scales. Scale identity
offsets inversely so the represented one-body correction remains unchanged.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from dfthcsos.verify import load_inputs,evaluate_words
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors,reconstruct_x
p=Path('initial_experiments/candidates/chemistry-c30-weighted.npz');z=dict(np.load(p));inp=load_inputs(p)
u=reconstruct_u_vectors(z['angle_words']/2048,z['tree_edges'])[76:].reshape(65,60,76)
X=reconstruct_x(u,inp.w);X=np.tril(X)+np.swapaxes(np.tril(X,k=-1),-1,-2)
i,j=np.tril_indices(76);multiplicity=np.where(i==j,1.,np.sqrt(2.))
Z=X[:,:,i,j].reshape(-1,len(i))*multiplicity
D=inp.ref['tei'][i[:,None],j[:,None],i[None,:],j[None,:]].copy()
D+=.5*inp.S[i,j,None]*(i==j)[None,:]
D+=.5*(i==j)[:,None]*inp.S[None,i,j]
D*=multiplicity[:,None]*multiplicity[None,:]
E=Z.T@Z-D
base_norm=float(np.linalg.norm(E));print('baseline tensor error',base_norm,flush=True)
a0=(Z*Z).sum(axis=1);unit=Z/np.sqrt(a0[:,None]);K=(unit@unit.T)**2;g=((unit@E)*unit).sum(axis=1)
print('computed normal equations',len(a0),flush=True)
base_metrics=evaluate_words(inp,z['angle_words'],z['tree_edges'],12);gap0=base_metrics['E_gap'];ww=inp.W.ravel()**2
calls=0

def objective(v):
 global calls
 t=1+v;delta=a0*v
 f=float(delta@(K@delta)+2*g@delta)
 grad=2*a0*(K@delta+g)
 gap=gap0+.5*np.sum(ww*(1/t-1))
 dg=-.5*ww/t**2
 penalty=1.
 if gap<3.92:
  f+=penalty*(gap-3.92)**2;grad+=2*penalty*(gap-3.92)*dg
 if gap>3.98:
  f+=penalty*(gap-3.98)**2;grad+=2*penalty*(gap-3.98)*dg
 calls+=1
 if calls%100==0:print(calls,'pred error',np.sqrt(max(0,base_norm**2+float(delta@(K@delta)+2*g@delta))),'gap',gap,'scale range',v.min(),v.max(),flush=True)
 return 1e6*f,1e6*grad


start_fit=dict(np.load('initial_experiments/candidates/chemistry-refit29.npz'))
start=(np.sum(start_fit['W']**2,axis=1)/np.sum(inp.w**2,axis=1)).ravel()-1
l1=np.abs(inp.w).sum(axis=1).ravel();signed=inp.w.sum(axis=1).ravel();offset=inp.W.ravel()
lambda1=base_metrics['Lambda']-.25*np.sum((l1+np.abs(offset+signed))**2)
calls=0

def effective_objective(v):
 global calls
 t=1+v;sqrt=np.sqrt(t);delta=a0*v;kd=K@delta
 err2=base_norm**2+delta@kd+2*g@delta
 derr=2*a0*(kd+g)
 gap=gap0+.5*np.sum(ww*(1/t-1));dgap=-.5*ww/t**2
 h=offset/sqrt+sqrt*signed
 A=sqrt*l1+np.abs(h)
 lam=lambda1+.25*np.sum(A*A)
 dlam=.5*A*(.5*l1/sqrt+np.sign(h)*(-.5*offset/(t*sqrt)+.5*signed/sqrt))
 rad=gap*(2*lam-gap)
 f=.5*np.log(rad);gradient=(gap*dlam+(lam-gap)*dgap)/rad
 limit=.0198**2
 if err2>limit:
  excess=err2/limit-1;f+=1000*excess**2;gradient+=2000*excess*derr/limit
 if gap<3.93:
  excess=(3.93-gap)/3.93;f+=1e5*excess**2;gradient-=2e5*excess*dgap/3.93
 calls+=1
 if calls%100==0:print('EFFECTIVE',calls,'error',np.sqrt(max(0,err2)),'gap',gap,'lambda_eff',np.sqrt(rad),flush=True)
 return f,gradient

res=minimize(effective_objective,start,jac=True,method='L-BFGS-B',bounds=[(-.15,.15)]*len(a0),options={'maxiter':900,'ftol':1e-14,'gtol':1e-8,'maxcor':20,'maxls':30})
print('OPT',res.message,flush=True)
t=(1+res.x).reshape(inp.W.shape);scale=np.sqrt(t);w=inp.w*scale[:,None,:];off=inp.W/scale
zz={**z,'W':w,'H2_I':off+w.sum(axis=1)}
np.savez('initial_experiments/candidates/chemistry-effective29.npz',**zz)
from dataclasses import replace
met=evaluate_words(replace(inp,params=zz,w=w,W=off),z['angle_words'],z['tree_edges'],12)
print('FINAL',json.dumps(met),flush=True)
