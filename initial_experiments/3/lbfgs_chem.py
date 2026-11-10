"""Polish a low-rank factorization with double-precision L-BFGS."""
import argparse,json,time
import numpy as np
import torch
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=300);p.add_argument('--columns',type=int)
p.add_argument('--float32',action='store_true');p.add_argument('--fixed-u',action='store_true')
p.add_argument('--adam',action='store_true');p.add_argument('--lr',type=float,default=2e-3)
a=p.parse_args();torch.set_num_threads(3);torch.backends.cuda.matmul.allow_tf32=False
dtype=torch.float32 if a.float32 else torch.float64;dev='cuda';HERE=Path(__file__).resolve().parent
z=np.load(a.init);u0=z['U'];w0=z['w'];W0=z['W'];R,B,C=w0.shape;n=76
ii,jj=np.tril_indices(n);f=np.where(ii==jj,1,np.sqrt(2.));ix=ii*n+jj
if a.columns and a.columns<C:
    pp=(u0[:,:,ii]*u0[:,:,jj])*f
    xx=np.einsum('rbc,rbp->rcp',w0,pp)
    _,vec=np.linalg.eigh(xx@xx.transpose(0,2,1));rot=vec[:,:,-a.columns:]
    w0=np.einsum('rbc,rcd->rbd',w0,rot);W0=np.einsum('rc,rcd->rd',W0,rot);C=a.columns
g=np.load(HERE/'parameters/chem_reference.npz')['G']
g=g[ix[:,None],ix[None,:]].astype(float)*f[:,None]*f[None,:]
G=torch.tensor(g,device=dev,dtype=dtype);fac=torch.tensor(f,device=dev,dtype=dtype)
ti=torch.tensor(ii,device=dev);tj=torch.tensor(jj,device=dev)
u=torch.tensor(u0,device=dev,dtype=dtype,requires_grad=not a.fixed_u)
w=torch.tensor(w0,device=dev,dtype=dtype,requires_grad=True)
def loss():
    un=u/u.norm(dim=2,keepdim=True)
    P=un[:,:,ti]*un[:,:,tj]*fac
    X=torch.bmm(w.transpose(1,2),P).reshape(-1,len(ii))
    D=X.T@X-G
    return D.square().sum()
params=[w] if a.fixed_u else [u,w]
opt=(torch.optim.Adam(params,lr=a.lr) if a.adam else
     torch.optim.LBFGS(params,lr=1,max_iter=10,history_size=40,line_search_fn='strong_wolfe',
                      tolerance_grad=1e-12,tolerance_change=1e-15))
best=float('inf');calls=0;start=time.time()
def closure():
    global calls
    opt.zero_grad();l=loss();l.backward();calls+=1
    return l
for s in range(a.steps):
    if a.adam:
        closure();opt.step()
        opt.param_groups[0]['lr']=a.lr*(.02+.98*.5*(1+np.cos(np.pi*s/a.steps)))
    else:opt.step(closure)
    with torch.no_grad():v=float(loss().sqrt())
    if v<best:
        best=v
        np.savez(a.out,U=(u/u.norm(dim=2,keepdim=True)).detach().cpu().numpy(),
                 w=w.detach().cpu().numpy(),W=W0,err2=v)
    if s%10==0 or s==a.steps-1:
        print(json.dumps(dict(step=s,calls=calls,error_2=v,best=best,seconds=time.time()-start)),flush=True)
