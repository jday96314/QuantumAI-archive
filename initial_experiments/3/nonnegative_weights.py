"""Convex nonnegative initialization against physical reference factors."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--reference',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=5000);a=p.parse_args();torch.set_num_threads(1)
z=dict(np.load(a.fit));ref=np.load(a.reference);u=z['U'][0];ur=ref['U'][0];wr=ref['w'][0];B=len(u)
ii,jj=np.tril_indices(76);fac=np.where(ii==jj,1,np.sqrt(2));P=u[:,ii]*u[:,jj]*fac
target=wr.T@(ur[:,ii]*ur[:,jj]*fac);d=dict(device='cuda',dtype=torch.float32)
pt=torch.tensor(P,**d);xt=torch.tensor(target,**d);S=pt@pt.T;H=pt@xt.T;step_size=1/torch.linalg.eigvalsh(S).max()
w=torch.zeros((B,len(target)),**d);y=w.clone();t=1.;began=time.time()
with torch.no_grad():
 for step in range(a.steps):
  new=(y-step_size*(S@y-H)).clamp(min=0);tn=.5*(1+np.sqrt(1+4*t*t));y=new+(t-1)/tn*(new-w);w=new;t=tn
  if step%500==0:print(json.dumps(dict(step=step,X_error=float((w.T@pt-xt).square().sum().sqrt()),seconds=time.time()-began)),flush=True)
z.pop('K',None);z.update(w=w.cpu().numpy().astype(float)[None,:,:],W=np.zeros((1,w.shape[1])))
np.savez(a.out,**z)
