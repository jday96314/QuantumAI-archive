"""ADMM feasibility acceleration between fixed native-factor Gram and PSD cone.

Fixed-factor polar projection preserves its tensor exactly. PSD projection
is generator-wise. The saved result records both tensor preservation and
error after PSD projection; it does not assert convergence or acceptance.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric
HERE=Path(__file__).resolve().parent
@torch.no_grad()
def main():
 p=argparse.ArgumentParser();p.add_argument('--init',type=Path,default=HERE/'parameters/psd_generator_rotation_native_C352.npz');p.add_argument('--steps',type=int,default=200);a=p.parse_args()
 torch.set_num_threads(1);start=time.time();z=dict(np.load(a.init));m=Metric(z['tree_edges']);X=m.tensor(z['factor']);Y=X.clone();dual=torch.zeros_like(X);rows=[]
 def project(T):
  ev,V=torch.linalg.eigh(m.unpack(T));M=(V*ev.clamp(min=0).unsqueeze(-2))@V.transpose(-1,-2);return M[:,m.ii,m.jj]*m.fac
 for step in range(a.steps+1):
  if step%10==0 or step==a.steps:
   positive=project(Y);ev=torch.linalg.eigvalsh(m.unpack(Y));rec=dict(step=step,negative_frobenius=float(ev.clamp(max=0).norm()),minimum_eigenvalue=float(ev.min()),positive_projection_error=float((positive.T@positive-m.G).norm()),seconds=time.time()-start)
   rows.append(rec);print(json.dumps(rec),flush=True)
  if step==a.steps:break
  Z=project(Y+dual);left,sv,right=torch.linalg.svd((Z-dual)@X.T,full_matrices=False);Y=(left@right)@X;dual+=Y-Z
 drift=float((Y.T@Y-X.T@X).norm());assert drift<1e-8,drift
 np.savez(HERE/'parameters/psd_generator_admm_native_C352.npz',factor=Y.cpu().numpy(),positive_factor=positive.cpu().numpy(),tree_edges=z['tree_edges'])
 out=dict(status='FINITE_PSD_GENERATOR_ADMM_SCREEN',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},rows=rows,tensor_drift=drift,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='fixed-tensor generator feasibility; no dictionary/native circuit')
 (HERE/'results/psd_generator_admm_native.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
