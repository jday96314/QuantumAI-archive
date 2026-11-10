"""Fixed exact-word chemistry rank reduction with physical gauge alignment."""
import os,sys,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for path in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:
 sys.path.insert(0,str(ROOT/path))
import numpy as np
import torch
from census_common import initial,quantities
torch.set_num_threads(2)
start=time.time();z,m,u,w,off,words=initial();P=m.feats(u)[0]
with torch.no_grad():
 Q,R=torch.linalg.qr(P.T);H=Q.T@m.G@Q;ev,V=torch.linalg.eigh((H+H.T)/2)
 X=w[0].T@P;outside=float((m.G-Q@H@Q.T).square().sum());norm=X.square().sum(1)
 rows=[]
 for C in [349,346,344,342,340]:
  wb=torch.linalg.solve(R,V[:,-C:]*ev[-C:].clamp(min=0).sqrt())
  Xb=wb.T@P
  for metric in ['physical_norm','weight_norm']:
   scores=norm if metric=='physical_norm' else w[0].square().sum(0)
   keep=torch.argsort(scores,descending=True)[:C].sort().values
   left,_,right=torch.linalg.svd(Xb@X[keep].T);rot=left@right
   wr=(wb@rot)[None];o=off[:,keep].clone();m.project_gap(u,wr,o)
   r=dict(C=C,removed_by=metric,**quantities(m,u,wr,o),elapsed=time.time()-start)
   rows.append(r);print(json.dumps(r),flush=True)
   Path('results/rank_trim_screen.json').write_text(json.dumps(rows,indent=2)+'\n')
   np.savez(f'parameters/rank_trim_C{C}_{metric}.npz',W=wr.cpu().numpy(),offset=o.cpu().numpy(),kept=keep.cpu().numpy(),source='banked exact words')
