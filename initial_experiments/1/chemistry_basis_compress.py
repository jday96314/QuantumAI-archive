"""Choose a smaller shared orbital-vector basis with weighted pivoted QR."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.linalg import qr
from dfthcsos.verify import load_inputs,evaluate_words,heuristic_costs
p=Path('from_host/sample_submission/dfthcsos/submission.npz');z=dict(np.load(p));inp=load_inputs(p)
# Preserve the decoded vectors exactly; regroup the angle-word rows with each basis.
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
u=reconstruct_u_vectors(z['angle_words']/2048,z['tree_edges'])[76:].reshape(65,60,76)
rows=np.tril_indices(76)
selected=[];projected=[]
for r in range(65):
 # Weighted rank-one matrices, with the lower triangle scaled to its Frobenius metric.
 mats=np.einsum('bp,bq->pqb',u[r],u[r])[rows].copy()
 mats[rows[0]!=rows[1]]*=np.sqrt(2)
 importance=np.linalg.norm(inp.w[r],axis=1)
 _,_,piv=qr(mats*importance[None,:],mode='economic',pivoting=True)
 selected.append(piv)
 # Solve for new coefficients against all original matrix columns.
 projected.append(mats)
for B in [58,56,52,48,40]:
 ww=[];hh=[];words=[z['angle_words'][:76]];uu=[]
 for r in range(65):
  ids=selected[r][:B];m=projected[r]
  # Include the offsets' effect on h1 by preserving its weight vector in the column space.
  transform=np.linalg.lstsq(m[:,ids],m,rcond=1e-12)[0]
  wr=transform@inp.w[r]
  ww.append(wr);hh.append(inp.W[r]+wr.sum(axis=0));uu.append(z['U'][r,ids]);words.append(z['angle_words'][76+r*60+ids])
 zz={**z,'U':np.array(uu),'W':np.array(ww),'H2_I':np.array(hh),'angle_words':np.concatenate(words)}
 out=Path(f'initial_experiments/candidates/chemistry-b{B}.npz');np.savez(out,**zz)
 candidate=replace(inp,params=zz,w=zz['W'],W=zz['H2_I']-zz['W'].sum(axis=1))
 metrics=evaluate_words(candidate,zz['angle_words'],z['tree_edges'],12)
 print('B',B,json.dumps(metrics),flush=True)
