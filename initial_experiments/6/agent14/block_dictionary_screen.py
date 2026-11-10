"""Screen exact-word R2 dictionaries via a block-covariance span lower bound.

The unrestricted (possibly indefinite, arbitrary-rank) optimum is a rigorous
numeric lower-bound diagnostic for these fixed words. It does not establish
physical validity and cannot exclude other dictionaries or repaired words.
"""
import sys,time,json
from pathlib import Path
ROOT=Path.cwd().parents[2]
for path in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:
 sys.path.insert(0,str(ROOT/path))
import numpy as np
import torch
from census_common import initial
start=time.time();z,m,u,w,off,words=initial();torch.set_num_threads(2)
with torch.no_grad():
 P=m.feats(u)[0];G=m.G;N=P.shape[0];K=w[0]@w[0].T
 # Spectral ordering of the incumbent covariance offers a distinct warm start
 # from training and then quantizing a new R2/B360 dictionary.
 graph=K.square();graph.fill_diagonal_(0);d=graph.sum(1).clamp(min=1e-30).rsqrt();graph=graph*d[:,None]*d[None,:]
 ev,V=torch.linalg.eigh(graph)
 orders={'original':torch.arange(N,device=P.device),'covariance_fiedler':torch.argsort(V[:,-2]),'physical_coordinate':torch.argsort(u[0,:,0])}
 gen=torch.Generator(device=P.device);gen.manual_seed(142030)
 for i in range(3):orders[f'random_{i}']=torch.randperm(N,generator=gen,device=P.device)
 rows=[]
 for name,order in orders.items():
  ids=[order[:N//2],order[N//2:]];Q=[torch.linalg.qr(P[ix].T)[0] for ix in ids]
  C=Q[0].T@Q[1];G1=Q[0].T@G@Q[0];G2=Q[1].T@G@Q[1]
  D=C@C.T;lam,V=torch.linalg.eigh(D);den=1-lam[:,None]*lam[None,:]
  A=G1-C@G2@C.T
  k1=V@((V.T@A@V)/den)@V.T;k2=G2-C.T@k1@C
  fitted=Q[0]@k1@Q[0].T+Q[1]@k2@Q[1].T
  err=float((G-fitted).norm());normal=max(float((Q[i].T@(G-fitted)@Q[i]).norm()) for i in range(2))
  r=dict(partition=name,block_sizes=[len(x) for x in ids],unrestricted_block_covariance_error=err,min_denominator=float(den.min()),normal_equation_residual=normal,below_20mHa=err<.02,elapsed=time.time()-start)
  rows.append(r);print(json.dumps(r),flush=True)
  Path('results/block_dictionary_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
  np.savez(f'parameters/block_partition_{name}.npz',order=order.cpu().numpy())
