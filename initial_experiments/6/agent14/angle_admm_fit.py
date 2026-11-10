"""Angle-lattice ADMM continuation with a native signed tensor objective.

The word constraint is optimized explicitly, rather than injected as fresh
noise. Saved integer fits remain screens until physical and native validation.
"""
import sys,argparse,json,time,os
from pathlib import Path
ROOT=Path.cwd().parents[2]
for pth in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.insert(0,str(ROOT/pth))
import numpy as np
import torch
from census_common import initial
from robust_fixed_fit import continuous_encode


def main():
 p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);p.add_argument('--outer',type=int,default=30);p.add_argument('--inner',type=int,default=30);p.add_argument('--rho',type=float,default=1e-8);p.add_argument('--growth',type=float,default=1.3);p.add_argument('--warm',type=int,default=200);p.add_argument('--angle-scale',type=float,default=1.);p.add_argument('--weight-ridge',type=float,default=1e-9);a=p.parse_args();torch.manual_seed(142105);torch.set_num_threads(2);start=time.time();z,m,_,_,_,_=initial();inp=dict(np.load(a.init));u=m.tensor(inp['U']);w=m.tensor(inp['w']);R,B,C=w.shape;coord=m.tensor(continuous_encode(u.cpu().numpy(),m.tree)/a.angle_scale).requires_grad_();eye=torch.eye(B,**m.d);P=m.feats(u);L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-10*eye);f=(L.transpose(1,2)@w).detach().requires_grad_();const=m.G.square().sum();dual=torch.zeros_like(coord);quant=(coord.detach()*a.angle_scale).round();rows=[];best=1e10
 def state():
  un=m.decode(coord*a.angle_scale);P=m.feats(un);L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-10*eye);wn=torch.linalg.solve_triangular(L.transpose(1,2),f,upper=True);X=(wn.transpose(1,2)@P).flatten(0,1);err=(const-2*(X*(X@m.G)).sum()+(X@X.T).square().sum()).clamp(min=1e-20);return un,wn,err
 def checkpoint(step,rho):
  nonlocal best
  with torch.no_grad():
   un,wn,err=state();X=(wn.transpose(1,2)@m.feats(un)).flatten(0,1);exact=float((X.T@X-m.G).norm());assert abs(exact**2-float(err))<1e-7;words=(coord*a.angle_scale).round();uq=m.decode(words);wq,qe=m.project(uq,wn,cycles=2);r=dict(step=step,rho=rho,continuous_error=exact,rounded_error=float(m.loss(uq,wn).sqrt()),refit_error=qe,rms_grid_distance=float((coord*a.angle_scale-(coord*a.angle_scale).round()).square().mean().sqrt()),rms_primal=float((coord*a.angle_scale-quant).square().mean().sqrt()),weight_norm=float(wn.norm()),seconds=time.time()-start);rows.append(r);print(json.dumps(r),flush=True)
   if qe<best:
    best=qe;np.savez(f'parameters/{a.tag}_best_integer.npz',U=uq.cpu().numpy(),w=wq.cpu().numpy(),angle_words=(words.cpu().numpy().astype(int)%256),tree_edges=m.tree,error2=qe)
   np.savez(f'parameters/{a.tag}_last.npz',U=un.cpu().numpy(),w=wn.cpu().numpy(),angles=(coord.detach()*a.angle_scale).cpu().numpy(),dual=dual.cpu().numpy(),quant=quant.cpu().numpy(),tree_edges=m.tree)
   Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',arguments=vars(a),rows=rows,finished=step==a.outer),indent=2)+'\n')
 checkpoint(-1,0)
 for outer in range(a.outer+1):
  rho=0 if outer==0 else a.rho*a.growth**(outer-1);opt=torch.optim.LBFGS([coord,f],lr=.7,max_iter=a.warm if outer==0 else a.inner,history_size=12,line_search_fn='strong_wolfe',tolerance_grad=1e-11,tolerance_change=1e-13)
  def closure():
   opt.zero_grad();un,wn,err=state();loss=err+rho*(coord*a.angle_scale-quant+dual).square().sum()+a.weight_ridge*wn.square().sum();loss.backward();assert torch.isfinite(coord.grad).all() and torch.isfinite(f.grad).all();return loss
  opt.step(closure)
  with torch.no_grad():
   quant=(coord*a.angle_scale+dual).round();dual.add_(coord.detach()*a.angle_scale-quant)
  checkpoint(outer,rho)
if __name__=='__main__':main()
