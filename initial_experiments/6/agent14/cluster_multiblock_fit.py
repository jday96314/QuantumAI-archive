"""Fit smaller physical blocks from localized incumbent covariance columns.

Native fixed target and exact word decoding remain unchanged. Results are
numerical screens until physical checks and complete native oracles pass.
"""
import sys,argparse,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for path in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.insert(0,str(ROOT/path))
import numpy as np
import torch
from scipy.optimize import linear_sum_assignment
from census_common import initial
from audit_core import encode


def initialize(m,u,w,R,B,C,varimax=False):
 W=w[0].clone();n=W.shape[1]
 if varimax:
  rotation=torch.eye(n,**m.d)
  for _ in range(40):
   V=W@rotation;left,_,right=torch.linalg.svd(W.T@(V.pow(3)-V*V.square().sum(0)[None,:]/len(V)),full_matrices=False);rotation=left@right
  W=W@rotation
 signatures=W.square().T;signatures=signatures/signatures.norm(dim=1,keepdim=True).clamp(min=1e-20);centers=[signatures[int(W.square().sum(0).argmax())]]
 while len(centers)<R:
  sim=signatures@torch.stack(centers).T;centers.append(signatures[int(sim.max(1).values.argmin())])
 centers=torch.stack(centers);capacity=(n+R-1)//R;labels=None
 for _ in range(30):
  score=signatures@centers.T;cost=-score[:,torch.arange(R,device=W.device).repeat_interleave(capacity)].cpu().numpy();ii,jj=linear_sum_assignment(cost);ll=np.empty(n,dtype=int);ll[ii]=jj//capacity
  if labels is not None and np.array_equal(labels,ll):break
  labels=ll;centers=torch.stack([signatures[torch.tensor(labels==r,device=W.device)].mean(0) for r in range(R)]);centers=centers/centers.norm(dim=1,keepdim=True)
 ids=[];weights=[]
 for r in range(R):
  cols=np.flatnonzero(labels==r);mass=W[:,cols].square().sum(1);ix=torch.argsort(mass,descending=True)[:B];ids.append(ix);ww=torch.zeros((B,C),**m.d);ww[:,:len(cols)]=W[ix][:,cols];weights.append(ww)
 un=torch.stack([u[0,ix] for ix in ids]);wn=torch.stack(weights);wn,error=m.project(un,wn,cycles=3);return un,wn,[ix.cpu().numpy() for ix in ids],error


def main():
 p=argparse.ArgumentParser();p.add_argument('--R',type=int,default=5);p.add_argument('--B',type=int,default=159);p.add_argument('--C',type=int,default=96);p.add_argument('--steps',type=int,default=2500);p.add_argument('--lr',type=float,default=.002);p.add_argument('--varimax',action='store_true');p.add_argument('--noise',type=int,default=0);p.add_argument('--ridge',type=float,default=1e-9);p.add_argument('--tag',default='cluster_R5_B159_C96');a=p.parse_args();start=time.time();torch.manual_seed(142083);torch.set_num_threads(2);z,m,u,w,off,words=initial();u,w,ids,initial_error=initialize(m,u,w,a.R,a.B,a.C,a.varimax);print(json.dumps(dict(stage='initialized',shape=list(w.shape),error=initial_error,seconds=time.time()-start)),flush=True);np.savez(f'parameters/{a.tag}_initial.npz',U=u.cpu().numpy(),w=w.cpu().numpy(),tree_edges=m.tree,ids=np.array(ids));coord=u.clone().requires_grad_();P=m.feats(u);eye=torch.eye(a.B,**m.d);L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-11*eye);f=(L.transpose(1,2)@w).detach().requires_grad_();opt=torch.optim.Adam([coord,f],lr=a.lr);const=m.G.square().sum();rows=[];best=1e10
 def state():
  un=coord/coord.norm(dim=-1,keepdim=True);P=m.feats(un);L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-11*eye);wn=torch.linalg.solve_triangular(L.transpose(1,2),f,upper=True);X=(wn.transpose(1,2)@P).flatten(0,1);err=(const-2*(X*(X@m.G)).sum()+(X@X.T).square().sum()).clamp(min=1e-20);return un,wn,X,err
 def noisy_state(un,perturbation):
  v=un;angles=torch.zeros((*un.shape[:-1],75),**m.d)
  for ts,ps,cs in reversed(m.stages):
   aa=v[...,ps];bb=v[...,cs];sign=torch.where(aa<0,-torch.ones_like(aa),torch.ones_like(aa));mag=aa.square()+bb.square();valid=mag>1e-24;safea=torch.where(valid,aa,torch.ones_like(aa));safeb=torch.where(valid,bb,torch.zeros_like(bb));angles[...,ts]=torch.atan2(safeb*sign,safea*sign)*(128/np.pi);new=v.clone();new[...,ps]=torch.where(valid,sign*mag.clamp(min=1e-24).sqrt(),torch.zeros_like(mag));new[...,cs]=0;v=new
  return m.decode(angles+perturbation)
 for step in range(a.steps+1):
  if step%250==0 or step==a.steps:
   with torch.no_grad():
    un,wn,X,err=state();exact=float((X.T@X-m.G).norm());assert abs(exact**2-float(err))<1e-7
    angle=encode(un.cpu().numpy().reshape(-1,76),m.tree,8,True);uq=m.decode(m.tensor(angle).reshape(a.R,a.B,75));wq,qe=m.project(uq,wn,cycles=2);cond=torch.linalg.cond(m.feats(uq)@m.feats(uq).transpose(1,2)).sqrt();r=dict(step=step,error=exact,rounded_error=float(m.loss(uq,wn).sqrt()),refit_error=qe,max_feature_condition=float(cond.max()),weight_norm=float(wn.norm()),seconds=time.time()-start);rows.append(r);print(json.dumps(r),flush=True)
    if exact<best:best=exact;np.savez(f'parameters/{a.tag}_best.npz',U=un.cpu().numpy(),w=wn.cpu().numpy(),tree_edges=m.tree,error2=exact)
    np.savez(f'parameters/{a.tag}_last_integer.npz',U=uq.cpu().numpy(),w=wq.cpu().numpy(),angle_words=angle,tree_edges=m.tree,error2=qe)
    Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',arguments=vars(a),finished=step==a.steps,rows=rows),indent=2)+'\n')
  if step==a.steps:break
  opt.zero_grad();un,wn,X,err=state();loss=err if not a.noise else err*0
  if a.noise:
   delta=torch.rand((*un.shape[:-1],75),**m.d)-.5
   for sample in range(a.noise):
    uq=noisy_state(un,delta if sample%2==0 else -delta);Pq=m.feats(uq);Xq=(wn.transpose(1,2)@Pq).flatten(0,1);loss=loss+(const-2*(Xq*(Xq@m.G)).sum()+(Xq@Xq.T).square().sum())/a.noise
  loss=loss+a.ridge*wn.square().sum();loss.backward();assert torch.isfinite(coord.grad).all() and torch.isfinite(f.grad).all();opt.step()
  with torch.no_grad():coord.div_(coord.norm(dim=-1,keepdim=True))
  opt.param_groups[0]['lr']=a.lr*(.05+.95*.5*(1+np.cos(np.pi*(step+1)/a.steps)))
if __name__=='__main__':main()
