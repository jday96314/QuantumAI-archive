"""Continuation from continuous Givens angles to an integer angle lattice."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=10);p.add_argument('--steps',type=int,default=12000)
p.add_argument('--lr',type=float,default=.02);p.add_argument('--ridge',type=float,default=1e-5)
p.add_argument('--columns',type=int);p.add_argument('--strength',type=float,default=100.)
p.add_argument('--start-strength',type=float,default=1e-4)
p.add_argument('--spectrum-rank',type=int);p.add_argument('--threads',type=int,default=3)
a=p.parse_args();torch.set_num_threads(a.threads);HERE=Path(__file__).resolve().parent;z=np.load(a.init)
u0=z['U'][0];B=len(u0);tree=z['tree_edges'] if 'tree_edges' in z else np.load(HERE/'candidates/chem35c45_rot.npz')['tree_edges']
v=u0.copy();th=np.zeros((B,75))
for j in reversed(range(75)):
 pa,ch=tree[j];th[:,j]=np.arctan2(v[:,ch],v[:,pa]);v[:,pa]=np.hypot(v[:,pa],v[:,ch]);v[:,ch]=0
alpha=np.pi/(1<<(a.brot-1));d=dict(device='cuda',dtype=torch.float64)
theta=torch.tensor(th/alpha,requires_grad=True,**d)
sp=np.load(HERE/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors']
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);const=G.square().sum()
spectral=None
if a.spectrum_rank:
 selected=np.argsort(ev)[-a.spectrum_rank:]
 spectral=torch.tensor(V[:,selected]*np.sqrt(np.maximum(ev[selected],0)),**d)
ti=torch.tensor(sp['ii'],device='cuda');tj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
eye=torch.eye(B,**d);startvec=torch.zeros((B,76),**d);startvec[:,0]=1
stages=[]
for j,(pa,ch) in enumerate(tree):
 pa,ch=int(pa),int(ch)
 if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
 stages[-1].append((j,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
def decode(t):
 u=startvec
 for ts,ps,cs in stages:
  ang=t[:,ts]*alpha;c,s=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
  new=u.clone();new[:,ps]=c*up-s*uc;new[:,cs]=s*up+c*uc;u=new
 return u
def state(t):
 u=decode(t);P=u[:,ti]*u[:,tj]*fac;S=(u@u.T).square()
 if spectral is None:H=P@G@P.T
 else:F=P@spectral;H=F@F.T
 Sr=S+a.ridge*eye
 with torch.no_grad():
  L=torch.linalg.cholesky(Sr)
  if a.columns and a.columns<B:
   hh=torch.linalg.solve_triangular(L,H,upper=False);hh=torch.linalg.solve_triangular(L,hh.T,upper=False).T
   e,v=torch.linalg.eigh(.5*(hh+hh.T));f=torch.linalg.solve_triangular(L.T,v[:,-a.columns:]*e[-a.columns:].clamp(min=0).sqrt(),upper=True);K=f@f.T
  else:K=torch.cholesky_solve(torch.cholesky_solve(H,L).T,L).T
  K=.5*(K+K.T)
 ks=K@Sr;loss=const-2*(K*H).sum()+(ks*ks.T).sum();ks=K@S
 error=(const-2*(K*H).sum()+(ks*ks.T).sum()).clamp(min=0).sqrt()
 return loss,error,u,K,S
with torch.no_grad():
 _,_,u,K,S=state(theta);yn=(K@S*K).sum(1).clamp(min=1e-12)
 # Exact coordinate derivative norms on the Givens tree, used as fixed penalty weights.
 cur=startvec.clone();J=torch.zeros((B,75,76),**d)
 for j,(pa,ch) in enumerate(tree):
  pa,ch=int(pa),int(ch);ang=theta[:,j]*alpha;c,s=ang.cos(),ang.sin()
  jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone();J[:,:,pa]=c[:,None]*jp-s[:,None]*jc;J[:,:,ch]=s[:,None]*jp+c[:,None]*jc
  up,uc=cur[:,pa].clone(),cur[:,ch].clone();cur[:,pa]=c*up-s*uc;cur[:,ch]=s*up+c*uc
  J[:,j,pa]=-alpha*cur[:,ch];J[:,j,ch]=alpha*cur[:,pa]
 weight=(8*J.square().sum(2)*yn[:,None]).clamp(min=1e-14)
opt=torch.optim.Adam([theta],lr=a.lr);best=float('inf');start=time.time()
for step in range(a.steps):
 opt.zero_grad();loss,error,u,K,S=state(theta)
 strength=10**(np.log10(a.start_strength)+(np.log10(a.strength)-np.log10(a.start_strength))*step/max(1,a.steps-1))
 penalty=(weight*(1-torch.cos(2*np.pi*theta))/(4*np.pi*np.pi)).sum()
 objective=loss+strength*penalty;objective.backward();opt.step()
 opt.param_groups[0]['lr']=a.lr*(.05+.95*.5*(1+np.cos(np.pi*step/a.steps)))
 if step%100==0 or step==a.steps-1:
  with torch.no_grad():
   words=theta.round();_,hard,uq,kq,sq=state(words);value=float(hard)
   if spectral is not None:
    pq=uq[:,ti]*uq[:,tj]*fac;hq=pq@G@pq.T;ks=kq@sq
    value=float((const-2*(kq*hq).sum()+(ks*ks.T).sum()).clamp(min=0).sqrt())
   if value<best:
    best=value;np.savez(a.out,U=uq.cpu().numpy()[None,:,:],K=kq.cpu().numpy(),angle_words=words.cpu().numpy().astype(np.int64)%(1<<a.brot),tree_edges=tree,brot=a.brot,err2=value)
   np.savez(str(Path(a.out).with_suffix(''))+'_soft.npz',U=u.detach().cpu().numpy()[None,:,:],theta=theta.cpu().numpy(),tree_edges=tree,err2=float(error))
  print(json.dumps(dict(step=step,soft_error=float(error.detach()),hard_error=value,best=best,strength=strength,penalty=float(penalty.detach()),seconds=time.time()-start)),flush=True)
