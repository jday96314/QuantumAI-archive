"""Spectral rank truncation of each physical generator into its own C1 block.

Unlike splitting rows, this preserves the original ensemble columns as warm
physical factors, while letting every factor have its own rank-B geometry.
"""
import sys,json,time,os
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
import numpy as np
import torch
from census_common import initial
from audit_core import save,encode
from dfthcsos.verify import heuristic_costs
start=time.time();z,m,u,w,off,words=initial();torch.set_num_threads(1)
with torch.no_grad():
 P=m.feats(u);F=(w.transpose(1,2)@P)[0];X=m.unpack(F);ev,V=torch.linalg.eigh(X);order=ev.abs().argsort(-1,descending=True);R=F.shape[0];rows=[]
 for B in [3,7,15,31,63,76]:
  ix=order[:,:B];weights=ev.gather(1,ix)[...,None];vectors=V.gather(2,ix[:,None,:].expand(-1,76,-1)).transpose(1,2);Xnew=(weights.transpose(1,2)@m.feats(vectors)).flatten(0,1);error=float((Xnew.T@Xnew-m.G).norm());aw=encode(vectors.cpu().numpy().reshape(-1,76),m.tree,8,True);uq=m.decode(m.tensor(aw).reshape(R,B,75));Xq=(weights.transpose(1,2)@m.feats(uq)).flatten(0,1);qe=float((Xq.T@Xq-m.G).norm());offset=off.reshape(R,1);g,l,s=m.physical(uq,weights,offset);hc=heuristic_costs(76,R,B,1,8);row=dict(R=R,B=B,C=1,float_error=error,integer_error=qe,Lambda=float(l),gap=float(g),lambda_eff=float(s.sqrt()),heuristic=hc,seconds=time.time()-start);rows.append(row);print(json.dumps(row),flush=True)
  if B in [7,15]:save(f'parameters/spectral_R{R}_B{B}_C1.npz',uq,weights,m.tree,aw,W=offset.cpu().numpy());save(f'parameters/spectral_R{R}_B{B}_C1_float.npz',vectors,weights,m.tree,W=offset.cpu().numpy())
 Path('results/rank_one_blocks_screen.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',rows=rows),indent=2)+'\n')
