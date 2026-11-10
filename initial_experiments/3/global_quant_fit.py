"""Quantization-aware global dictionary fitting with exact linear solves."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=12);p.add_argument('--steps',type=int,default=8000)
p.add_argument('--lr',type=float,default=.03);p.add_argument('--ridge',type=float,default=1e-5)
p.add_argument('--psd',action='store_true');p.add_argument('--columns',type=int)
p.add_argument('--preconditioned',action='store_true');p.add_argument('--jitter',type=float,default=0.)
a=p.parse_args();torch.set_num_threads(3);HERE=Path(__file__).resolve().parent
z=np.load(a.init);u0=z['U'][0];B=len(u0);tree=z['tree_edges'] if 'tree_edges' in z else np.load(HERE/'candidates/chem35c45_rot.npz')['tree_edges']
sp=np.load(HERE/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors'];ii,jj=sp['ii'],sp['jj'];fac=sp['fac']
d=dict(device='cuda',dtype=torch.float64);G=torch.tensor((V*(np.maximum(ev,0) if a.psd else ev))@V.T,**d);const=G.square().sum()
theta=torch.tensor(encode(u0,tree,a.brot),**d,requires_grad=True)
if a.jitter:
    torch.manual_seed(20261023)
    with torch.no_grad():theta.add_((torch.rand_like(theta)*2-1)*a.jitter)
ti=torch.tensor(ii,device='cuda');tj=torch.tensor(jj,device='cuda');ff=torch.tensor(fac,**d)
stages=[]
for t,(pa,ch) in enumerate(tree):
    pa,ch=int(pa),int(ch)
    if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
    stages[-1].append((t,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
startvec=torch.zeros((B,76),**d);startvec[:,0]=1
def decode():
    rounded=theta.round();angles=(theta+(rounded-theta).detach())*(np.pi/(1<<(a.brot-1)))
    u=startvec
    for ts,ps,cs in stages:
        ang=angles[:,ts];co,sn=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
        new=u.clone();new[:,ps]=co*up-sn*uc;new[:,cs]=sn*up+co*uc;u=new
    return u,rounded
eye=torch.eye(B,**d);opt=(torch.optim.SGD([theta],lr=a.lr,momentum=.9) if a.preconditioned else torch.optim.Adam([theta],lr=a.lr));best=1e9;start=time.time();hdiag=None
def raw_curvature(K,S):
    v=startvec.clone();J=torch.zeros((B,75,76),**d);alpha=np.pi/(1<<(a.brot-1))
    for t,(pa,ch) in enumerate(tree):
        pa,ch=int(pa),int(ch);ang=theta[:,t].round()*alpha;co,sn=ang.cos(),ang.sin()
        jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone();J[:,:,pa]=co[:,None]*jp-sn[:,None]*jc;J[:,:,ch]=sn[:,None]*jp+co[:,None]*jc
        vp,vc=v[:,pa].clone(),v[:,ch].clone();v[:,pa]=co*vp-sn*vc;v[:,ch]=sn*vp+co*vc
        J[:,t,pa]=-alpha*v[:,ch];J[:,t,ch]=alpha*v[:,pa]
    yn=(K@S*K).sum(1).clamp(min=1e-14)
    return (8*J.square().sum(2)*yn[:,None]).clamp(min=1e-14)
for step in range(a.steps):
    opt.zero_grad();u,words=decode();P=u[:,ti]*u[:,tj]*ff
    S=(u@u.T).square();H=(P@G)@P.T;Sr=S+a.ridge*eye
    with torch.no_grad():
        L=torch.linalg.cholesky(Sr)
        if a.columns and a.columns<B:
            hh=torch.linalg.solve_triangular(L,H,upper=False);hh=torch.linalg.solve_triangular(L,hh.T,upper=False).T
            ee,vv=torch.linalg.eigh(.5*(hh+hh.T));f=torch.linalg.solve_triangular(L.T,vv[:,-a.columns:]*ee[-a.columns:].clamp(min=0).sqrt(),upper=True);K=f@f.T
        else:K=torch.cholesky_solve(torch.cholesky_solve(H,L).T,L).T
        K=.5*(K+K.T)
    KS=K@Sr;loss=const-2*(K*H).sum()+(KS*KS.T).sum();loss.backward()
    if a.preconditioned:
        with torch.no_grad():
            if hdiag is None or step%20==0:hdiag=raw_curvature(K,S)
            theta.grad.div_(hdiag).clamp_(-.25,.25)
    opt.step()
    with torch.no_grad():
        physical=K@S;error=float((const-2*(K*H).sum()+(physical*physical.T).sum()).clamp(min=0).sqrt())
    if step%5==0 and error<best:
        best=error
        np.savez(a.out,U=u.detach().cpu().numpy()[None,:,:],K=K.cpu().numpy(),err2=error,
                 angle_words=words.detach().cpu().numpy().astype(np.int64)%(1<<a.brot),tree_edges=tree,brot=a.brot)
    if step%100==0:print(json.dumps(dict(step=step,error_2=error,best=best,lambda_gram_bound=float(.25*K.abs().sum()),seconds=time.time()-start)),flush=True)
    opt.param_groups[0]['lr']=a.lr*(.1+.9*.5*(1+np.cos(np.pi*step/a.steps)))
z=dict(np.load(a.out));e,v=np.linalg.eigh(z.pop('K'));order=np.argsort(e)[::-1];e=np.maximum(e[order],0);v=v[:,order]
z.update(w=(v*np.sqrt(e))[None,:,:],W=np.zeros((1,B)));np.savez(a.out,**z)
