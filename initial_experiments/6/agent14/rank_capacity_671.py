"""Exact rank/normalization screens for the newly repaired integer B671 words."""
import sys,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
z,u,w,tree,words=load('parameters/lll_B671_native.npz');m=Metric(tree);torch.set_num_threads(2);ut=m.tensor(u);old=m.tensor(w);off=m.tensor(z['W']);start=time.time()
with torch.no_grad():
 P=m.feats(ut)[0];Q,T=torch.linalg.qr(P.T);H=Q.T@m.G@Q;ev,V=torch.linalg.eigh((H+H.T)/2);outside=(m.G-Q@H@Q.T).square().sum();rows=[]
 for C in [350,360,380,400,450,500,671]:
  wb=torch.linalg.solve(T,V[:,-C:]*ev[-C:].clamp(min=0).sqrt());warm=torch.zeros((671,C),**m.d);warm[:,:min(C,350)]=old[0,:,:min(C,350)];l,_,r=torch.linalg.svd(wb.T@warm,full_matrices=False);wb=wb@(l@r);ww=wb[None];oo=torch.zeros((1,C),**m.d);oo[:,:min(C,350)]=off[:,:min(C,350)];m.project_gap(ut,ww,oo);g,la,s=m.physical(ut,ww,oo);er=float((outside+(H-(V[:,-C:]*ev[-C:].clamp(min=0))@V[:,-C:].T).square().sum()).sqrt());row=dict(C=C,error=er,Lambda=float(la),gap=float(g),lambda_eff=float(s.sqrt()),seconds=time.time()-start);rows.append(row);print(json.dumps(row),flush=True)
  if er<.02:save(f'parameters/rank_671_C{C}.npz',ut,ww,tree,words,W=oo.cpu().numpy())
 Path('results/rank_capacity_671.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',rows=rows),indent=2)+'\n')
