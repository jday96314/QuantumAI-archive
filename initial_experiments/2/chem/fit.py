"""GPU refit of a DF-THC factorization at reduced (R, B, C): stage 1 minimizes error_2."""
import argparse,time,numpy as np,torch
p=argparse.ArgumentParser()
p.add_argument('--ref',required=True);p.add_argument('--init',default=None)
p.add_argument('--R',type=int);p.add_argument('--B',type=int);p.add_argument('--C',type=int)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--lr',type=float,default=3e-3)
p.add_argument('--out',required=True);p.add_argument('--lam',type=float,default=0.0);p.add_argument('--seed',type=int,default=0)
a=p.parse_args()
torch.backends.cuda.matmul.allow_tf32=False;dev='cuda'
ref=np.load(a.ref);G=torch.tensor(ref['G'],device=dev)
if a.init:
    z=np.load(a.init);U0=z['U'];w0=z['w'];W0=z['W']
else:
    U0=ref['U'];w0=ref['w'];W0=ref['W']
R0,B0,C0=w0.shape
# Prune: keep strongest rho, then strongest b per rho, then strongest c.
rs=np.argsort(-np.abs(w0).sum(axis=(1,2)))[:a.R]
U=[];w=[];W=[]
for r in rs:
    bs=np.argsort(-np.abs(w0[r]).sum(axis=1))[:a.B]
    U.append(U0[r,bs]);w.append(w0[r][bs]);W.append(W0[r])
U=np.array(U);w=np.array(w);W=np.array(W)
if a.C<w.shape[2]:
    cs=np.argsort(-np.abs(w).sum(axis=(0,1)))[:a.C];w=w[:,:,cs];W=W[:,cs]
elif a.C>w.shape[2]:
    rng=np.random.default_rng(a.seed);extra=a.C-w.shape[2]
    w=np.concatenate([w,1e-3*rng.standard_normal(w.shape[:2]+(extra,))],2);W=np.concatenate([W,np.zeros((W.shape[0],extra))],1)
u=torch.tensor(U,device=dev,dtype=torch.float32,requires_grad=True)
wt=torch.tensor(w,device=dev,dtype=torch.float32,requires_grad=True)
def err2():
    un=u/u.norm(dim=2,keepdim=True)
    X=torch.einsum('rbp,rbq,rbc->rcpq',un,un,wt).reshape(-1,76*76)
    D=X.T@X-G
    return D.norm()
Wfix=torch.tensor(W,device=dev,dtype=torch.float32)
def lam2():
    wB=Wfix+wt.sum(1)
    return 0.25*((wt.abs().sum(1)+wB.abs())**2).sum()
opt=torch.optim.Adam([u,wt],lr=a.lr)
sched=torch.optim.lr_scheduler.CosineAnnealingLR(opt,a.steps)
t0=time.time();best=1e9
for s in range(a.steps):
    opt.zero_grad();e=err2();(e*e+a.lam*lam2()).backward();opt.step();sched.step()
    v=e.item()
    if v<best:
        best=v;un=(u/u.norm(dim=2,keepdim=True)).detach().cpu().numpy().astype(np.float64)
        np.savez(a.out,U=un,w=wt.detach().cpu().numpy().astype(np.float64),W=W,err2=v)
    if s%200==0 or s==a.steps-1:print(s,'err2 %.5f best %.5f'%(v,best),'%.0fs'%(time.time()-t0),flush=True)
