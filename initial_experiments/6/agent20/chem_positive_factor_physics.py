"""Intrinsic physical offset screen for a PSD generator factor.

For positive dictionaries, the coefficient L1 norm equals each matrix trace.
Thus Lambda/gap can be computed directly from the factor, independently of
any dictionary construction. Tensor accuracy alone is not native acceptance.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric
HERE=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--init',type=Path,default=HERE/'parameters/psd_generator_admm_native_C352.npz');a=p.parse_args();torch.set_num_threads(1);start=time.time();z=dict(np.load(a.init));m=Metric(z['tree_edges'],device='cpu');X=m.tensor(z['positive_factor']);trace=X[:,m.ii==m.jj].sum(1);o=(-trace).detach().requires_grad_()
 def physical(off):
  h=m.h0-m.unpack(off@X);e=torch.linalg.eigvalsh(h);ab=e.abs().sum();g=m.constant+ab-e.sum()+.5*off.square().sum();l=ab+.25*(trace+(off+trace).abs()).square().sum();return g,l,g*(2*l-g)
 def record(name):
  with torch.no_grad():
   g,l,e=physical(o);return dict(name=name,gap=float(g),Lambda=float(l),lambda_eff=float(e.sqrt()) if e>0 else None,offset_norm=float(o.norm()))
 rows=[record('negative_trace')];opt=torch.optim.LBFGS([o],lr=1,max_iter=400,tolerance_grad=1e-10,tolerance_change=1e-13,line_search_fn='strong_wolfe')
 def closure():
  opt.zero_grad();g=physical(o)[0];g.backward();return g
 opt.step(closure);rows.append(record('minimum_gap'))
 if float(physical(o)[0])<3.91515:
  low=o.detach().clone();high=torch.zeros_like(o)
  assert float(physical(high)[0])>3.91515
  for _ in range(80):
   mid=(low+high)/2
   if float(physical(mid)[0])<3.91515:low=mid
   else:high=mid
  with torch.no_grad():o.copy_(high)
  rows.append(record('floor_gap_interpolation'))
 with torch.no_grad():
  ev=torch.linalg.eigvalsh(m.unpack(X));err=float((X.T@X-m.G).norm())
 np.savez(HERE/'parameters/positive_factor_physics.npz',factor=X.cpu().numpy(),W=o.detach().numpy(),tree_edges=z['tree_edges'])
 out=dict(status='INTRINSIC_PSD_FACTOR_PHYSICS',error2=err,minimum_generator_eigenvalue=float(ev.min()),rows=rows,seconds=time.time()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='positive-factor trace normalization under native formulas; no finite dictionary or oracle')
 (HERE/'results/positive_factor_physics.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
