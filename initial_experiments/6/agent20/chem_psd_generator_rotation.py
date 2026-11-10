"""Finite PSD-generator feasibility screen at a fixed native tensor factor.

Alternates spectral projection of each generator with exact orthogonal
Procrustes rotation of a fixed tensor square root. Rotation preserves X.T X;
PSD cone projection is used only as a target. Zero-padding tests additional
generator counts without changing the source tensor.
"""
import argparse,time,json,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric
from chem_rank_screen import BANK
HERE=Path(__file__).resolve().parent

@torch.no_grad()
def main():
 p=argparse.ArgumentParser();p.add_argument('--counts',type=int,nargs='+',default=[350,700]);p.add_argument('--steps',type=int,default=200);p.add_argument('--native-rank',type=int,default=0);a=p.parse_args()
 torch.set_num_threads(1);started=time.time();z=dict(np.load(BANK));m=Metric(z['tree_edges']);u=m.decode(m.tensor(z['angle_words'][76:]))
 X=m.tensor(z['W'])[0].T@m.feats(u)
 if a.native_rank:
  spectrum,eigenvectors=torch.linalg.eigh(m.G);X=(eigenvectors[:,-a.native_rank:]*spectrum[-a.native_rank:].clamp(min=0).sqrt()).T
 K=X.shape[0];factor_error=float((X.T@X-m.G).norm());reports=[]
 # Trace vector rotation makes all output traces equal and positive.
 diag=m.ii==m.jj;t=X[:,diag].sum(1)
 for C in a.counts:
  if C<K:raise ValueError('count must retain the source rank')
  trace=torch.zeros(C,**m.d);trace[:K]=t;v=trace/trace.norm()-torch.ones(C,**m.d)/C**.5
  H=torch.eye(C,**m.d)-2*torch.outer(v,v)/v.square().sum();A=H[:,:K];Y=A@X;history=[]
  for step in range(a.steps+1):
   matrices=m.unpack(Y);ev,V=torch.linalg.eigh(matrices);neg=ev.clamp(max=0)
   if step%10==0 or step==a.steps:
    rec=dict(count=C,step=step,negative_frobenius=float(neg.norm()),minimum_eigenvalue=float(ev.min()),psd_generators=int((ev[:,0]>=-1e-10).sum()),seconds=time.time()-started)
    history.append(rec);print(json.dumps(rec),flush=True)
   if step==a.steps:break
   pos=(V*ev.clamp(min=0).unsqueeze(-2))@V.transpose(-1,-2)
   target=pos[:,m.ii,m.jj]*m.fac
   # Rectangular polar factor is an isometry even for C>K.
   left,sv,right=torch.linalg.svd(target@X.T,full_matrices=False);A=left@right;Y=A@X
  pos=(V*ev.clamp(min=0).unsqueeze(-2))@V.transpose(-1,-2)
  positive=pos[:,m.ii,m.jj]*m.fac;positive_error=float((positive.T@positive-m.G).norm())
  tensor_drift=float((Y.T@Y-X.T@X).norm());assert tensor_drift<1e-8,tensor_drift
  np.savez(HERE/'parameters'/f'psd_generator_rotation_{"native" if a.native_rank else "incumbent"}_C{C}.npz',factor=Y.cpu().numpy(),rotation=A.cpu().numpy(),positive_factor=positive.cpu().numpy(),tree_edges=z['tree_edges'])
  reports.append(dict(C=C,history=history,tensor_drift=tensor_drift,positive_projection_error=positive_error))
  out=dict(status='FINITE_FIXED_TENSOR_PSD_GENERATOR_ROTATION',arguments=vars(a),source_factor_error=factor_error,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reports=reports,seconds=time.time()-started,scope='finite feasibility search; no dictionary or native circuit')
  (HERE/f'results/psd_generator_rotation_{"native" if a.native_rank else "incumbent"}.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
