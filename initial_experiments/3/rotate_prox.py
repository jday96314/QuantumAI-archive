"""Sparse proximal maps alternated with exact orthogonal Procrustes projection."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--gamma',type=float,default=.01)
a=p.parse_args();torch.set_num_threads(3);z=dict(np.load(a.fit));w=z['w'];W=z['W'];R,B,C=w.shape;assert R==1
A=np.concatenate([w[0],(W[0]+w[0].sum(0))[None,:]],0)
at=torch.tensor(A,device='cuda',dtype=torch.float32);q=torch.eye(C,device='cuda');cur=at.clone()
best=float(.25*np.square(abs(A).sum(0)).sum());bestq=np.eye(C);start=time.time()
ks=torch.arange(1,len(A)+1,device='cuda')[:,None]
with torch.no_grad():
 for step in range(a.steps):
  gamma=a.gamma*(.001+.999*(1-step/a.steps)**2)
  d=cur.abs();ds=d.sort(dim=0,descending=True).values
  thresholds=2*gamma*ds.cumsum(0)/(1+2*gamma*ks)
  active=(ds>thresholds).sum(0).clamp(min=1)-1
  tau=thresholds.gather(0,active[None,:])[0]
  y=cur.sign()*(d-tau).clamp(min=0)
  u,s,vh=torch.linalg.svd(at.T@y,full_matrices=False);q=u@vh;cur=at@q
  value=float(.25*cur.abs().sum(0).square().sum())
  if value<best:best=value;bestq=q.cpu().numpy().astype(float)
  if step%100==0:print(json.dumps(dict(step=step,lambda_h2=value,best=best,gamma=gamma,seconds=time.time()-start)),flush=True)
q,_=np.linalg.qr(bestq);out=A@q
assert np.allclose(A@A.T,out@out.T,rtol=1e-10,atol=1e-10)
z['w']=out[:B][None,:,:];z['W']=(out[B]-out[:B].sum(0))[None,:]
np.savez(a.out,**z);print('final',float(.25*np.square(abs(out).sum(0)).sum()),flush=True)
