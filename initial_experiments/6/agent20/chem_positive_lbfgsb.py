"""Bound-constrained native tensor fit of positive rank-one dictionaries.

Direct nonnegative coefficients admit activation at zero. L-BFGS-B uses
float64 exact native residual/gradients instead of log weights or float32
Adam. Rounding is evaluated separately, never differentiated through.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import minimize
from audit_core import Metric,load,save,encode
HERE=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True);p.add_argument('--iterations',type=int,default=500);p.add_argument('--tag',required=True);a=p.parse_args();torch.set_num_threads(1);start=time.time();z,u,w,tree,_=load(a.init);m=Metric(tree);R,B,C=w.shape;nu=u.size;shapeu=u.shape;shapew=w.shape;scale=.1
 x=np.r_[u.ravel(),w.ravel()/scale];bounds=[(None,None)]*nu+[(0,None)]*w.size;rows=[];calls=0;last={};best=float('inf')
 def unpack(v):
  u=v[:nu].reshape(shapeu);return u/u.norm(dim=-1,keepdim=True),scale*v[nu:].reshape(shapew)
 def fun(x):
  nonlocal calls,last
  v=m.tensor(x).requires_grad_();uu,ww=unpack(v);X=(ww.transpose(1,2)@m.feats(uu)).flatten(0,1);loss=(X.T@X-m.G).square().sum()*1000
  loss.backward();value=float(loss.detach());grad=v.grad.detach().cpu().numpy().copy();calls+=1;last=dict(error2=(value/1000)**.5,calls=calls,seconds=time.time()-start)
  return value,grad
 def callback(x):
  nonlocal best
  if len(rows)%25==0:print(json.dumps(dict(iteration=len(rows)+1,**last)),flush=True)
  rows.append(dict(iteration=len(rows)+1,**last))
  if last['error2']<best:
   best=last['error2'];uu,ww=unpack(m.tensor(x));save(HERE/'parameters'/f'{a.tag}_float.npz',uu,ww,tree,W=(-ww.sum(1)).cpu().numpy())
 result=minimize(fun,x,jac=True,method='L-BFGS-B',bounds=bounds,callback=callback,options=dict(maxiter=a.iterations,maxls=40,maxcor=10,ftol=1e-15,gtol=1e-8))
 with torch.no_grad():
  uu,ww=unpack(m.tensor(result.x));words=encode(uu.cpu().numpy().reshape(-1,76),tree,8,True);uq=m.decode(m.tensor(words).reshape(R,B,75))
  err=float(m.loss(uu,ww).sqrt());rounded=float(m.loss(uq,ww).sqrt());save(HERE/'parameters'/f'{a.tag}_rounded.npz',uq,ww,tree,words,W=(-ww.sum(1)).cpu().numpy())
 out=dict(status='FINITE_POSITIVE_LBFGSB',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},optimizer_success=bool(result.success),termination=str(result.message),rows=rows,final_error2=err,rounded_error2=rounded,calls=calls,seconds=time.time()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='native float64 tensor optimization with nonnegative coefficients; rounded checkpoint is not an accepted submission')
 (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='rows'}),flush=True)
if __name__=='__main__':main()
