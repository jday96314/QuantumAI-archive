"""Shared positive rank-one dictionaries for the native PSD generator factor.

Warm directions are common-generator mean eigenvectors and spectral atoms
of individual generators. Nonnegative coefficients use NNLS multiplicative
updates before smooth joint vector/weight fitting to the generator matrices.
Float and exact 8-bit dictionaries are independently evaluated against the
full native tensor. This is a new representation screen, no circuit claim.
"""
import argparse,copy,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric,encode,save
HERE=Path(__file__).resolve().parent

@torch.no_grad()
def nnls(P,target,w,steps):
 S=P@P.transpose(-1,-2);T=P@target.transpose(-1,-2)
 for _ in range(steps):w*=T.clamp(min=1e-20)/(S@w).clamp(min=1e-20)
 return w

def main():
 p=argparse.ArgumentParser();p.add_argument('--init',type=Path,default=HERE/'parameters/psd_generator_admm_native_C352.npz');p.add_argument('--R',type=int,default=4);p.add_argument('--B',type=int,default=240);p.add_argument('--steps',type=int,default=4000);p.add_argument('--tag',default='positive_R4B240');p.add_argument('--checkpoint',type=int,default=500);p.add_argument('--lr',type=float,default=.005);p.add_argument('--tensor-loss',action='store_true');p.add_argument('--init-dictionary',type=Path);p.add_argument('--direct-weights',action='store_true');a=p.parse_args()
 torch.set_num_threads(1);torch.manual_seed(2026102524);started=time.time();z=dict(np.load(a.init));m=Metric(z['tree_edges']);X=m.tensor(z['positive_factor']);C=len(X)//a.R;assert len(X)%a.R==0;target=X.reshape(a.R,C,-1)
 with torch.no_grad():
  matrices=m.unpack(target);_,common=torch.linalg.eigh(matrices.mean(1));common=common.transpose(-1,-2)
  ev,V=torch.linalg.eigh(matrices);strength=ev[:,:,-1];pool=V[:,:,:,-4:].transpose(-1,-2).reshape(a.R,-1,76)
  order=ev[:,:,-4:].reshape(a.R,-1).argsort(-1,descending=True);extra=max(0,a.B-76)
  # Weighted systematic sampling covers each generator, rather than only the
  # nearly identical principal eigendirections of the common positive mean.
  chosen=torch.linspace(0,pool.shape[1]-1,extra,device=X.device).round().long()
  u=torch.cat([common,pool[:,chosen]],1) if extra else common[:,-a.B:]
  w=torch.ones((a.R,a.B,C),**m.d)*target[:, :, m.ii==m.jj].sum(-1).mean()/a.B
  w=nnls(m.feats(u),target,w,1000)
 if a.init_dictionary:
  dz=dict(np.load(a.init_dictionary));u=m.tensor(dz['U']);w=m.tensor(dz['w']);assert w.shape==(a.R,a.B,C)
 # Smooth generator-coordinate objective; no rounding estimator.
 mt=copy.copy(m);mt.d=dict(device='cuda',dtype=torch.float32);mt.fac=m.fac.float();mt.G=m.G.float()
 uu=u.float().detach().requires_grad_(True);logw=(w.float() if a.direct_weights else w.float().clamp(min=1e-8).log()).detach().requires_grad_(True);tg=target.float()
 opt=torch.optim.Adam([uu,logw],lr=a.lr);rows=[];best=float('inf')
 for step in range(a.steps+1):
  if step%a.checkpoint==0 or step==a.steps:
   with torch.no_grad():
    un=uu.double();un/=un.norm(dim=-1,keepdim=True);ww=logw.double().clamp(min=0) if a.direct_weights else logw.double().exp();P=m.feats(un);pred=ww.transpose(1,2)@P
    error=float((pred.flatten(0,1).T@pred.flatten(0,1)-m.G).norm());projection_error=float((pred-target).norm())
    words=encode(un.cpu().numpy().reshape(-1,76),z['tree_edges'],8,True);uq=m.decode(m.tensor(words).reshape(a.R,a.B,75));wq=nnls(m.feats(uq),pred if a.tensor_loss else target,ww.clone(),1000)
    qp=wq.transpose(1,2)@m.feats(uq);qe=float((qp.flatten(0,1).T@qp.flatten(0,1)-m.G).norm())
    rec=dict(step=step,shape=[a.R,a.B,C],float_error=error,exact_word_error=qe,factor_error=projection_error,min_weight=float(ww.min()),weight_norm=float(ww.norm()),seconds=time.time()-started)
    rows.append(rec);print(json.dumps(rec),flush=True)
    save(HERE/'parameters'/f'{a.tag}_float.npz',un,ww,z['tree_edges'],W=(-ww.sum(1)).cpu().numpy())
    save(HERE/'parameters'/f'{a.tag}_integer.npz',uq,wq,z['tree_edges'],words,W=(-wq.sum(1)).cpu().numpy())
    if qe<best:
     best=qe;save(HERE/'parameters'/f'{a.tag}_best_integer.npz',uq,wq,z['tree_edges'],words,W=(-wq.sum(1)).cpu().numpy())
    out=dict(status='FINITE_POSITIVE_DICTIONARY_FIT',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},rows=rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),finished=step==a.steps,scope='native target error with nonnegative generator coefficients; no physical offset optimization or native circuit')
    (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(out,indent=2)+'\n')
  if step==a.steps:break
  opt.zero_grad();un=uu/uu.norm(dim=-1,keepdim=True);ww=logw.clamp(min=0) if a.direct_weights else logw.exp();P=mt.feats(un);pred=ww.transpose(1,2)@P
  if a.tensor_loss:
   pf=pred.flatten(0,1);loss=(pf.T@pf-mt.G).square().sum()
  else:loss=(pred-tg).square().sum()
  loss.backward();opt.step()
  with torch.no_grad():uu/=uu.norm(dim=-1,keepdim=True);logw.clamp_(0,100) if a.direct_weights else logw.clamp_(-23,8)
  lr=a.lr*(.01**((step+1)/a.steps))
  for g in opt.param_groups:g['lr']=lr
if __name__=='__main__':main()
