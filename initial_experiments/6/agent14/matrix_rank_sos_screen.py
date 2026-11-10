"""Independent low matrix-rank SOS blocks from incumbent tensor factors.

Screen only: each factor X_c gets its own small orthogonal rank-one dictionary.
No accuracy, normalization, or compiled-oracle claim is made for a screen.
"""
import sys,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for path in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.insert(0,str(ROOT/path))
import numpy as np
import torch
from census_common import initial
from audit_core import encode
from dfthcsos.verify import heuristic_costs


def main():
 start=time.time();z,m,u,w,off,words=initial();rows=[]
 with torch.no_grad():
  X=(w.transpose(1,2)@m.feats(u)).flatten(0,1);mat=m.unpack(X);ev,V=torch.linalg.eigh(mat);order=torch.argsort(ev.abs(),descending=True);ev=torch.gather(ev,1,order);V=torch.gather(V,2,order[:,None,:].expand(-1,76,-1));R=len(X)
  for B in [1,2,3,4,6,8,12,16,24,32,48,76]:
   uu=V[:,:,:B].transpose(1,2);ww=ev[:,:B,None];err=float(m.loss(uu,ww).sqrt());h=heuristic_costs(76,R,B,1,8);r=dict(R=R,B=B,C=1,float_error=err,heuristic_T=h['toffoli'],heuristic_Q=h['qubits']+32,angle_rows=76+R*B,seconds=time.time()-start)
   if B<=8:
    angles=encode(uu.cpu().numpy().reshape(-1,76),m.tree,8,True);uq=m.decode(m.tensor(angles).reshape(R,B,75));r['rounded_error']=float(m.loss(uq,ww).sqrt());np.savez(f'parameters/matrix_rank_B{B}_initial.npz',U=uu.cpu().numpy(),w=ww.cpu().numpy(),tree_edges=m.tree,error2=err)
   rows.append(r);print(json.dumps(r),flush=True)
  Path('results/matrix_rank_sos_screen.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',rows=rows),indent=2)+'\n')
if __name__=='__main__':main()
