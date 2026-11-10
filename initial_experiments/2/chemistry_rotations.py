"""Optimize ensemble column coordinates without changing the represented tensors."""
import json
from pathlib import Path
import numpy as np
from dfthcsos.verify import load_inputs,evaluate_words
from dfthcsos.scoreinfra.dfthc_math import lambda_h2
p=Path('from_host/sample_submission/dfthcsos/submission.npz')
z=dict(np.load(p));w=z['W'];h=z['H2_I'];W=h-w.sum(axis=1)
print('baseline lambda2',lambda_h2(w,W),flush=True)
ww=np.empty_like(w);hh=np.empty_like(h)
for r in range(w.shape[0]):
 u,s,v=np.linalg.svd(np.vstack([w[r],h[r]]),full_matrices=False)
 ww[r]=w[r]@v.T;hh[r]=h[r]@v.T
print('svd lambda2',lambda_h2(ww,hh-ww.sum(axis=1)),flush=True)
print('singular example',s,flush=True)
for r in [0,10,30,60]:
 s=np.linalg.svd(np.vstack([w[r],h[r]]),compute_uv=False)
 print('rank',r, np.round(s,7).tolist(),flush=True)
# Pairwise rotations minimize column l1-norm sum of squares, independently for rho.
# Include H2_I in the l1 norm: lambda2=1/4 sum ||[w;h] column||_1^2.
a=np.concatenate([w,h[:,None,:]],axis=1).copy()
rng=np.random.default_rng(20261023)
for sweep in range(20):
 improvement=0
 for i,j in rng.permutation(np.array([(i,j) for i in range(60) for j in range(i)])):
  x=a[:,:,i].copy();y=a[:,:,j].copy()
  before=(np.abs(x).sum(axis=1)**2+np.abs(y).sum(axis=1)**2)
  best=before.copy();best_angle=np.zeros(65)
  # Coarse coordinate descent; fine refinement will follow promising runs.
  for t in np.linspace(-np.pi/4,np.pi/4,17):
   c,s=np.cos(t),np.sin(t)
   val=np.abs(c*x+s*y).sum(axis=1)**2+np.abs(-s*x+c*y).sum(axis=1)**2
   take=val<best;best[take]=val[take];best_angle[take]=t
  c=np.cos(best_angle)[:,None];s=np.sin(best_angle)[:,None]
  a[:,:,i]=c*x+s*y;a[:,:,j]=-s*x+c*y
  improvement+=float((before-best).sum()/4)
 print('sweep',sweep,'lambda2',float(np.sum(np.abs(a).sum(axis=1)**2)/4),'improvement',improvement,flush=True)
 np.savez('initial_experiments/candidates/chemistry-rotated.npz',**{**z,'W':a[:,:60,:],'H2_I':a[:,60,:]})
