"""Measure ensemble SVD truncation against the host's unchanged tensor metrics."""
import json
from pathlib import Path
import numpy as np
from dfthcsos.verify import load_inputs,evaluate_words,heuristic_costs
p=Path('from_host/sample_submission/dfthcsos/submission.npz');z=dict(np.load(p));inp=load_inputs(p)
a=np.concatenate([z['W'],z['H2_I'][:,None,:]],axis=1).astype(float)
u,s,v=np.linalg.svd(a,full_matrices=False)
print('row weight quantiles',np.quantile(np.linalg.norm(z['W'],axis=2),[0,.01,.1,.5,.9,1]),flush=True)
for C in [40,30,25,20,15]:
 aa=u[:,:,:C]*s[:,None,:C]
 zz={**z,'W':aa[:,:60,:],'H2_I':aa[:,60,:]}
 path=Path(f'initial_experiments/candidates/chemistry-c{C}.npz');np.savez(path,**zz)
 candidate=load_inputs(zz)
 met=evaluate_words(candidate,z['angle_words'],z['tree_edges'],12)
 print('C',C,json.dumps(met),flush=True)
