"""Rotate fixed physical factors toward a nonnegative weight factorization."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--threshold',type=float,default=0.)
a=p.parse_args();torch.set_num_threads(2);z=dict(np.load(a.fit));w=z['w'];W=z['W'];R,B,C=w.shape;assert R==1
A=w[0].copy();wb=W[0]+A.sum(0);at=torch.tensor(A,device='cuda',dtype=torch.float32)
q=torch.eye(C,device='cuda');bestq=np.eye(C);best=float(.25*np.square(abs(A).sum(0)+abs(wb)).sum());began=time.time()
with torch.no_grad():
 for step in range(a.steps):
  rotated=at@q;sign=torch.where(rotated.sum(0)<0,-1.,1.);q=q*sign;rotated=rotated*sign
  target=(rotated-a.threshold*(1-step/a.steps)).clamp(min=0)
  left,_,right=torch.linalg.svd(at.T@target,full_matrices=False);q=left@right;rotated=at@q
  shifted=torch.tensor(wb,device='cuda',dtype=torch.float32)@q
  value=float(.25*(rotated.abs().sum(0)+shifted.abs()).square().sum())
  if value<best:best=value;bestq=q.cpu().numpy().astype(float)
  if step%100==0:print(json.dumps(dict(step=step,lambda_h2=value,best=best,negative_squared=float(rotated.clamp(max=0).square().sum()),seconds=time.time()-began)),flush=True)
q,_=np.linalg.qr(bestq);out=A@q;off=wb@q-out.sum(0)
assert np.allclose(A@A.T,out@out.T,rtol=1e-10,atol=1e-10)
z['w']=out[None,:,:];z['W']=off[None,:];np.savez(a.out,**z)
print(json.dumps(dict(final_lambda_h2=float(.25*np.square(abs(out).sum(0)+abs(wb@q)).sum()))),flush=True)
