"""Rotate a fixed-word rank-compressed factor to reduce its L1 normalization."""
import sys,json,time,subprocess
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.insert(0,str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
import numpy as np
import torch
from census_common import initial
from jacobi_l1 import encode,decode
from audit_core import save
z,m,u,original,off,words=initial();torch.set_num_threads(2)
for C in [340]:
 r=np.load(f'parameters/rank_trim_C{C}_physical_norm.npz');w=r['W'];offset=r['offset'];tag=f'rank_gauge_C{C}';rows=[]
 for sweep in range(4):
  A=np.concatenate([w[0],(offset[0]+w[0].sum(0))[None]],axis=0);src=Path(f'parameters/{tag}_input.bin');dst=src.with_name(f'{tag}_output.bin');encode(A,src)
  with Path(f'results/{tag}_{sweep}.log').open('w') as log:subprocess.run([str(ROOT/'initial_experiments/4/codex/jacobi_l1'),str(src),str(dst),'10','12000',str(142126+sweep)],stdout=log,check=True)
  B=decode(dst);assert np.max(np.abs(A@A.T-B@B.T))<1e-9;w=B[None,:-1];offset=(B[-1]-B[:-1].sum(0))[None];wt=m.tensor(w);ot=m.tensor(offset);m.project_gap(u,wt,ot);offset=ot.cpu().numpy();g,l,s=m.physical(u,wt,ot);row=dict(sweep=sweep,error=float(m.loss(u,wt).sqrt()),gap=float(g),Lambda=float(l),lambda_eff=float(s.sqrt()));rows.append(row);print(json.dumps(row),flush=True);save(f'parameters/{tag}.npz',u,wt,m.tree,words,W=offset);Path(f'results/{tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',rows=rows),indent=2)+'\n')
