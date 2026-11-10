"""Optimize identity offsets under the strict SOS-gap floor, preserving h2."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from dfthcsos.verify import load_inputs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors,reconstruct_x,build_h1_prime
p=Path('initial_experiments/candidates/chemistry-c36-refined.npz');z=dict(np.load(p));inp=load_inputs(p)
uall=reconstruct_u_vectors(z['angle_words']/2048,z['tree_edges']);u=uall[76:].reshape(65,60,76)
X=reconstruct_x(u,inp.w)
X=np.tril(X)+np.swapaxes(np.tril(X,k=-1),-1,-2)
flat=X.reshape(-1,76*76)
h1=inp.ref['oei']-.5*np.einsum('prrq->pq',inp.ref['tei'])
base=h1-.5*113*inp.S+inp.beta_1*np.eye(76)
weights=inp.w.sum(axis=1).ravel();l1=np.abs(inp.w).sum(axis=1).ravel()
start=inp.W.ravel()
floor=3.9155
counter=0

def objective(w,penalty):
 global counter
 h=base-(w@flat).reshape(76,76)
 ev,V=np.linalg.eigh(h);lam1=np.abs(ev).sum()
 hh=w+weights
 lam2=.25*np.sum((l1+np.abs(hh))**2)
 lam=lam1+lam2
 gap=float(inp.energy_upper-inp.ref['ecore']+113*inp.beta_1+(lam1-ev.sum())+.5*(w@w))
 signmatrix=(V*np.sign(ev))@V.T
 gradlam= -flat@signmatrix.ravel()+.5*(l1+np.abs(hh))*np.sign(hh)
 gradgap= -flat@(signmatrix-np.eye(76)).ravel()+w
 # Smooth objective for the effective normalization, squared to avoid roots.
 f=gap*(2*lam-gap)
 g=2*(gap*gradlam+(lam-gap)*gradgap)
 if gap<floor:
  f+=penalty*(floor-gap)**2
  g-=2*penalty*(floor-gap)*gradgap
 counter+=1
 if counter%20==0:print(counter,'lambda',lam,'gap',gap,'lambda_eff',np.sqrt(max(0,gap*(2*lam-gap))),flush=True)
 return f,g

w=start.copy()
for penalty in [1e4,1e5,1e6,1e7]:
 res=minimize(lambda v:objective(v,penalty),w,jac=True,method='L-BFGS-B',options={'maxiter':150,'ftol':1e-12,'gtol':1e-7,'maxls':30,'maxcor':12})
 w=res.x
 print('OPT',penalty,res.message,flush=True)
 out={**z,'H2_I':w.reshape(inp.W.shape)+inp.w.sum(axis=1)}
 np.savez('initial_experiments/candidates/chemistry-offsets.npz',**out)
