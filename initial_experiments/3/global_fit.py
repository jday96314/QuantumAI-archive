"""Variable-projection fit of one global THC dictionary.

At each vector update, solve the complete optimal linear coefficient matrix
with a small ridge. The envelope derivative avoids differentiating the solve.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=5000);p.add_argument('--lr',type=float,default=.002)
p.add_argument('--ridge',type=float,default=1e-7)
p.add_argument('--optimizer',choices=['adam','lbfgs'],default='adam')
p.add_argument('--columns',type=int,help='constrain the projected coefficient matrix rank')
p.add_argument('--psd',action='store_true',help='fit the positive spectral part of the reference')
p.add_argument('--spectrum-rank',type=int,help='use a spectral approximation for gradients; checkpoints use the full target')
a=p.parse_args();torch.set_num_threads(3);HERE=Path(__file__).resolve().parent
z=np.load(a.init);u0=z['U'][0];B=len(u0);sp=np.load(HERE/'parameters/reference_spectrum.npz')
ev,V=sp['eigenvalues'],sp['eigenvectors'];ii,jj=sp['ii'],sp['jj'];fac=sp['fac']
d=dict(device='cuda',dtype=torch.float64)
G=torch.tensor((V*(np.maximum(ev,0) if a.psd else ev))@V.T,**d);const=G.square().sum()
full_const=const
Z=None
if a.spectrum_rank:
    assert a.psd
    ne=min(a.spectrum_rank,len(ev));ix=np.argsort(ev)[-ne:]
    Z=torch.tensor(V[:,ix]*np.sqrt(np.maximum(ev[ix],0)),**d)
    const=torch.tensor(np.square(np.maximum(ev[ix],0)).sum(),**d)
u=torch.tensor(u0,**d,requires_grad=True);ti=torch.tensor(ii,device='cuda');tj=torch.tensor(jj,device='cuda');ff=torch.tensor(fac,**d)
eye=torch.eye(B,**d)
opt=(torch.optim.Adam([u],lr=a.lr) if a.optimizer=='adam' else
     torch.optim.LBFGS([u],lr=a.lr,max_iter=12,history_size=30,line_search_fn='strong_wolfe',
                       tolerance_grad=1e-12,tolerance_change=1e-15))
best=1e9;start=time.time();current={}
def closure():
    opt.zero_grad();un=u/u.norm(dim=1,keepdim=True)
    P=un[:,ti]*un[:,tj]*ff
    S=(un@un.T).square()
    if Z is None:H=(P@G)@P.T
    else:F=P@Z;H=F@F.T
    Sr=S+a.ridge*eye
    with torch.no_grad():
        L=torch.linalg.cholesky(Sr)
        if a.columns and a.columns<B:
            Vh=torch.linalg.solve_triangular(L,H,upper=False)
            Vh=torch.linalg.solve_triangular(L,Vh.T,upper=False).T
            ee,vv=torch.linalg.eigh(.5*(Vh+Vh.T))
            factor=vv[:,-a.columns:]*ee[-a.columns:].clamp(min=0).sqrt()
            factor=torch.linalg.solve_triangular(L.T,factor,upper=True)
            K=factor@factor.T
        else:K=torch.cholesky_solve(torch.cholesky_solve(H,L).T,L).T
        K=.5*(K+K.T)
    KS=K@Sr
    loss=const-2*(K*H).sum()+(KS*KS.T).sum()
    loss.backward()
    with torch.no_grad():
        physical=K@S
        error=float((const-2*(K*H).sum()+(physical*physical.T).sum()).clamp(min=0).sqrt())
    current.update(un=un.detach(),K=K,error=error)
    return loss
for step in range(a.steps):
    if a.optimizer=='adam':closure();opt.step()
    else:opt.step(closure)
    un,K,error=current['un'],current['K'],current['error']
    if Z is not None and step%5==0:
        with torch.no_grad():
            P=un[:,ti]*un[:,tj]*ff;S=(un@un.T).square();H=(P@G)@P.T;ks=K@S
            error=float((full_const-2*(K*H).sum()+(ks*ks.T).sum()).clamp(min=0).sqrt())
    if step%10==0 and error<best:
        best=error
        np.savez(a.out,U=un.detach().cpu().numpy()[None,:,:],K=K.cpu().numpy(),err2=error)
    if step%(100 if a.optimizer=='adam' else 5)==0:
        print(json.dumps(dict(step=step,error_2=error,best=best,
                              lambda_gram_bound=float(.25*K.abs().sum()),seconds=time.time()-start)),flush=True)
    if a.optimizer=='adam':opt.param_groups[0]['lr']=a.lr*(.1+.9*.5*(1+np.cos(np.pi*step/a.steps)))
z=dict(np.load(a.out));e,v=np.linalg.eigh(z.pop('K'));order=np.argsort(e)[::-1];e=np.maximum(e[order],0);v=v[:,order]
C=a.columns or B
z.update(w=(v[:,:C]*np.sqrt(e[:C]))[None,:,:],W=np.zeros((1,C)))
np.savez(a.out,**z)
