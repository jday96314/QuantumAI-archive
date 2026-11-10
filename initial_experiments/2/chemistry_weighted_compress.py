"""Compress ensemble columns using the decoded matrix norm, not coefficient norm."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from dfthcsos.verify import load_inputs,evaluate_words
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=Path('from_host/sample_submission/dfthcsos/submission.npz');z=dict(np.load(p));inp=load_inputs(p)
u=reconstruct_u_vectors(z['angle_words']/2048,z['tree_edges'])[76:].reshape(65,60,76)
G=(u@u.transpose(0,2,1))**2
for gamma in [.1,1,10]:
 gram=inp.w.transpose(0,2,1)@G@inp.w+gamma*inp.W[:,:,None]*inp.W[:,None,:]
 ev,vec=np.linalg.eigh(gram);vec=vec[:,:,::-1]
 for C in [24,28,32,34]:
  rotation=vec[:,:,:C]
  ww=inp.w@rotation;off=np.einsum('rc,rck->rk',inp.W,rotation)
  zz={**z,'W':ww,'H2_I':off+ww.sum(axis=1)}
  candidate=replace(inp,params=zz,w=ww,W=off)
  met=evaluate_words(candidate,z['angle_words'],z['tree_edges'],12)
  name=f'chemistry-weighted-g{gamma}-c{C}'
  np.savez(f'initial_experiments/candidates/{name}.npz',**zz)
  print(name,json.dumps(met),flush=True)
