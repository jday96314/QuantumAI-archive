"""Optimize the complete column rotation with differentiable QR retraction.

Optimization uses fp32; the final orthogonal matrix is retracted in fp64 and
applied to the original fp64 factors, preserving their tensor accurately.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=4000);p.add_argument('--lr',type=float,default=.001)
p.add_argument('--seed',type=int);p.add_argument('--mix',type=float,default=1.)
p.add_argument('--smooth',type=float,default=0.)
p.add_argument('--mix-penalty',type=float,default=0.,help='sign-mixing coefficient relative to the H2 normalization objective')
a=p.parse_args();torch.set_num_threads(3);torch.backends.cuda.matmul.allow_tf32=False
z=dict(np.load(a.fit));w=z['w'];W=z['W'];R,B,C=w.shape;assert R==1
A=np.concatenate([w[0],(W[0]+w[0].sum(0))[None,:]],0)
at=torch.tensor(A,device='cuda',dtype=torch.float32)
o=torch.eye(C,device='cuda')
if a.seed is not None:
    torch.manual_seed(a.seed);o=o+a.mix*torch.randn_like(o)/np.sqrt(C)
o.requires_grad_(True)
opt=torch.optim.Adam([o],lr=a.lr);best=.25*np.square(abs(A).sum(0)).sum();bestq=np.eye(C);start=time.time()
for s in range(a.steps):
    opt.zero_grad();q,_=torch.linalg.qr(o);rot=at@q
    eps=a.smooth*(1-s/a.steps)**2
    loss=.25*(torch.sqrt(rot.square()+eps*eps) if eps else rot.abs()).sum(0).square().sum()
    mixing=(rot[:B].abs().sum(0).square()-rot[:B].sum(0).square()).sum()
    if a.mix_penalty:loss=loss+a.mix_penalty*mixing
    value=float((.25*rot.abs().sum(0).square().sum()).detach())
    if value<best:best=value;bestq=q.detach().cpu().numpy().astype(float)
    loss.backward();opt.step();opt.param_groups[0]['lr']=a.lr*(.02+.98*.5*(1+np.cos(np.pi*s/a.steps)))
    if s%100==0:print(json.dumps(dict(step=s,lambda_h2=value,best=best,mixing=float(mixing.detach()),mix_penalty=a.mix_penalty,seconds=time.time()-start)),flush=True)
q,_=np.linalg.qr(bestq);out=A@q
assert np.allclose(A@A.T,out@out.T,rtol=1e-10,atol=1e-10)
z['w']=out[:B][None,:,:];z['W']=(out[B]-out[:B].sum(0))[None,:]
np.savez(a.out,**z);print('final',float(.25*np.square(abs(out).sum(0)).sum()),flush=True)
