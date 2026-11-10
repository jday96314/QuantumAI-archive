"""Refine the orthogonal ensemble coordinates using fine Jacobi rotations."""
from pathlib import Path
import numpy as np
p=Path('initial_experiments/candidates/chemistry-rotated.npz');z=dict(np.load(p))
a=np.concatenate([z['W'],z['H2_I'][:,None,:]],axis=1).astype(float)
R,D,C=a.shape;rng=np.random.default_rng(20261024)
for sweep in range(100):
 scale=[.10,.04,.015,.006,.002,.0008,.0003,.0001][min(sweep//12,7)]
 improvement=0
 for i,j in rng.permutation(np.array([(i,j) for i in range(C) for j in range(i)])):
  x=a[:,:,i].copy();y=a[:,:,j].copy()
  before=(np.abs(x).sum(axis=1)**2+np.abs(y).sum(axis=1)**2)
  best=before.copy();best_angle=np.zeros(R)
  for t in np.linspace(-scale,scale,13):
   c,s=np.cos(t),np.sin(t)
   val=np.abs(c*x+s*y).sum(axis=1)**2+np.abs(-s*x+c*y).sum(axis=1)**2
   take=val<best;best[take]=val[take];best_angle[take]=t
  c=np.cos(best_angle)[:,None];s=np.sin(best_angle)[:,None]
  a[:,:,i]=c*x+s*y;a[:,:,j]=-s*x+c*y
  improvement+=float((before-best).sum()/4)
 print('sweep',sweep,'lambda2',float(np.sum(np.abs(a).sum(axis=1)**2)/4),'improvement',improvement,flush=True)
 np.savez('initial_experiments/candidates/chemistry-refined.npz',**{**z,'W':a[:,:D-1,:],'H2_I':a[:,D-1,:]})
