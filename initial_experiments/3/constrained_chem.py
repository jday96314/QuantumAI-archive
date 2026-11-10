"""Reduce the actual chemistry normalization with an augmented Lagrangian."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=12000);p.add_argument('--cap',type=float,default=.0175)
p.add_argument('--lr',type=float,default=1e-4)
a=p.parse_args();torch.set_num_threads(3);torch.backends.cuda.matmul.allow_tf32=False
HERE=Path(__file__).resolve().parent;dev='cuda';z=np.load(a.source)
if 'w' in z:
    w0=z['w'];off0=z['W'];U=z['U']
else:
    w0=z['W'];off0=z['H2_I']-w0.sum(1)
    R,B,C=w0.shape
    U=reconstruct_u_vectors(z['angle_words'][76:]/2048,z['tree_edges']).reshape(R,B,76)
R,B,C=w0.shape;n=76
ref=np.load(HERE/'parameters/chem_reference.npz');ii,jj=np.tril_indices(n);f=np.where(ii==jj,1,np.sqrt(2.));ix=ii*n+jj
G=ref['G'][ix[:,None],ix[None,:]].astype(float)*f[:,None]*f[None,:]
G=torch.tensor(G,device=dev,dtype=torch.float32)
u=torch.tensor(U,device=dev,dtype=torch.float32,requires_grad=True)
w=torch.tensor(w0,device=dev,dtype=torch.float32,requires_grad=True)
W=torch.tensor(off0,device=dev,dtype=torch.float64,requires_grad=True)
fac=torch.tensor(f,device=dev,dtype=torch.float32)
ti=torch.tensor(ii,device=dev);tj=torch.tensor(jj,device=dev)
h0=torch.tensor(ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(n),device=dev,dtype=torch.float64)
gapconst=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
def values():
    un=u/u.norm(dim=2,keepdim=True)
    P=un[:,:,ti]*un[:,:,tj]*fac
    X=torch.bmm(w.transpose(1,2),P).reshape(-1,len(ii))
    e2=(X.T@X-G).square().sum()
    un64=un.double();ww=w.double()
    v=(ww*W[:,None,:]).sum(2)
    h=h0-torch.einsum('rb,rbp,rbq->pq',v,un64,un64)
    eig=torch.linalg.eigvalsh(.5*(h+h.T));lh1=eig.abs().sum()
    eps=1e-6
    aw=torch.sqrt(ww.square()+eps*eps).sum(1)
    wb=W+ww.sum(1)
    lh2=.25*(aw+torch.sqrt(wb.square()+eps*eps)).square().sum()
    lam=lh1+lh2
    gap=gapconst+lh1-eig.sum()+.5*W.square().sum()
    le2=gap*(2*lam-gap)
    return e2,lam,gap,le2
opt=torch.optim.Adam([{'params':[u,w],'lr':a.lr},{'params':[W],'lr':a.lr*5}])
mue=0.;mug=0.;rho=2.;best=1e9;start=time.time()
for s in range(a.steps):
    opt.zero_grad();e2,lam,gap,le2=values()
    ge=e2/(a.cap*a.cap)-1
    gg=(3.916-gap)/3.916
    def pen(g,mu):return (torch.relu(mu+rho*g).square()-mu*mu)/(2*rho)
    loss=le2/2500+pen(ge,mue)+pen(gg,mug)
    loss.backward();opt.step()
    if s%50==0:
        # Score and arrays refer to the same post-step state.
        with torch.no_grad():ee,ll,gp,l2=values()
        error=float(ee.sqrt());lv=float(l2.clamp(min=0).sqrt())
        if error<=a.cap*1.005 and float(gp)>=3.915 and lv<best:
            best=lv
            np.savez(a.out,U=(u/u.norm(dim=2,keepdim=True)).detach().cpu().numpy().astype(float),
                     w=w.detach().cpu().numpy().astype(float),W=W.detach().cpu().numpy(),err2=error,lambda_eff=lv)
        if s%100==0:
            mue=max(0.,mue+rho*float(ge));mug=max(0.,mug+rho*float(gg))
            print(json.dumps(dict(step=s,error_2=error,Lambda=float(ll),gap=float(gp),lambda_eff=lv,
                                  best=best,mue=mue,mug=mug,seconds=time.time()-start)),flush=True)
    lr=a.lr*(.1+.9*(1-s/a.steps))
    opt.param_groups[0]['lr']=lr;opt.param_groups[1]['lr']=lr*5
