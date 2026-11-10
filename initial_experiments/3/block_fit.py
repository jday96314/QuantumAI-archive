"""Fit compact multi-ensemble THC shapes with a stable symmetric feature loss."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=30000);p.add_argument('--lr',type=float,default=.001)
p.add_argument('--ridge',type=float,default=1e-7);p.add_argument('--normalize-trace',action='store_true');a=p.parse_args();torch.set_num_threads(3)
torch.backends.cuda.matmul.allow_tf32=False;HERE=Path(__file__).resolve().parent
z=np.load(a.init);u0,w0=z['U'],z['w'];sp=np.load(HERE/'parameters/reference_spectrum.npz')
d=dict(device='cuda',dtype=torch.float32);ev,V=sp['eigenvalues'],sp['eigenvectors']
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda')
fac=torch.tensor(sp['fac'],**d);u=torch.tensor(u0,requires_grad=True,**d);w=torch.tensor(w0,requires_grad=True,**d)
opt=torch.optim.Adam([u,w],lr=a.lr);best=1e9;start=time.time()
iv=(sp['ii']==sp['jj']).astype(float);vt=V.T@iv;trace_norm=float(np.maximum(ev,0)@(vt*vt))
for step in range(a.steps):
 opt.zero_grad();un=u/u.norm(dim=2,keepdim=True);P=un[:,:,ii]*un[:,:,jj]*fac
 ww=w*torch.sqrt(trace_norm/w.sum(1).square().sum().clamp(min=1e-15)) if a.normalize_trace else w
 X=torch.bmm(ww.transpose(1,2),P).reshape(-1,len(sp['ii']));res=X.T@X-G;e2=res.square().sum()
 reg=ww.square().sum();loss=e2+a.ridge*reg;loss.backward()
 err=float(e2.detach().sqrt())
 if step%50==0 and err<best:
  best=err;np.savez(a.out,U=un.detach().cpu().numpy().astype(float),w=ww.detach().cpu().numpy().astype(float),
                          W=np.zeros((w.shape[0],w.shape[2])),err2=err)
 if step%200==0:print(json.dumps(dict(step=step,error_2=err,best=best,weight_norm=float(reg.detach()),seconds=time.time()-start)),flush=True)
 opt.step()
 opt.param_groups[0]['lr']=a.lr*(.03+.97*.5*(1+np.cos(np.pi*step/a.steps)))
