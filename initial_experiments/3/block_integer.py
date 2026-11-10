"""Integer Givens descent alternating with PSD block variable projection."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode

p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--cycles',type=int,default=12)
p.add_argument('--moves',type=int,default=16);p.add_argument('--brot',type=int,default=8)
a=p.parse_args();torch.set_num_threads(2);torch.backends.cuda.matmul.allow_tf32=False
root=Path(__file__).resolve().parent;z=dict(np.load(a.init));R,B,N=z['U'].shape;C=z['w'].shape[2]
tree=z.get('tree_edges',np.load(root/'candidates/chem35c45_rot.npz')['tree_edges'])
sp=np.load(root/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors']
d=dict(device='cuda',dtype=torch.float64);alpha=np.pi/(1<<(a.brot-1));T=len(tree)
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);const=G.square().sum()
ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
words=z.get('angle_words')
if words is None or int(z.get('brot',a.brot))!=a.brot:words=encode(z['U'].reshape(-1,N),tree,a.brot)
theta=torch.tensor(np.asarray(words).reshape(R*B,T),requires_grad=True,**d)
w=torch.tensor(z['w'],**d);start=torch.zeros((R*B,N),**d);start[:,0]=1
stages=[]
for j,(pa,ch) in enumerate(tree):
 pa,ch=int(pa),int(ch)
 if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
 stages[-1].append((j,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
def decode():
 u=start
 for ts,ps,cs in stages:
  ang=theta[:,ts]*alpha;c,s=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs]
  new=u.clone();new[:,ps]=c*up-s*uc;new[:,cs]=s*up+c*uc;u=new
 return u.reshape(R,B,N)
def features(u):return u[:,:,ii]*u[:,:,jj]*fac
@torch.no_grad()
def project(u,w,cycles):
 P=features(u);S=P@P.transpose(1,2);L=torch.linalg.cholesky(S+1e-11*torch.eye(B,**d))
 Q=torch.linalg.solve_triangular(L,P,upper=False);flat=Q.reshape(R*B,-1)
 ov=(flat@flat.T).reshape(R,B,R,B).permute(0,2,1,3).contiguous()
 cross=(flat@G@flat.T).reshape(R,B,R,B);target=torch.stack([cross[r,:,r,:] for r in range(R)])
 f=L.transpose(1,2)@w;H=f@f.transpose(1,2)
 residual=target-torch.einsum('rsab,sbc,rsdc->rad',ov,H,ov)
 for _ in range(cycles):
  for r in range(R):
   h=.5*(residual[r]+residual[r].T)+H[r];ee,vv=torch.linalg.eigh(h)
   ff=vv[:,-C:]*ee[-C:].clamp(min=0).sqrt();new=ff@ff.T;change=new-H[r];H[r]=new
   A=ov[:,r];residual-=A@change@A.transpose(1,2)
 ee,vv=torch.linalg.eigh(H);ff=vv[:,:,-C:]*ee[:,-C:].clamp(min=0).sqrt()[:,None,:]
 return torch.linalg.solve_triangular(L.transpose(1,2),ff,upper=True)
def loss(u,w):
 P=features(u);X=(w.transpose(1,2)@P).reshape(R*C,-1)
 # Feature-space contractions avoid a full 2926-square residual for every trial.
 return const-2*((X@G)*X).sum()+(X@X.T).square().sum()
@torch.no_grad()
def curvature(u,w):
 v=start.clone();J=torch.zeros((R*B,T,N),**d)
 for j,(pa,ch) in enumerate(tree):
  pa,ch=int(pa),int(ch);ang=theta[:,j]*alpha;c,s=ang.cos(),ang.sin()
  jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone();J[:,:,pa]=c[:,None]*jp-s[:,None]*jc;J[:,:,ch]=s[:,None]*jp+c[:,None]*jc
  up,uc=v[:,pa].clone(),v[:,ch].clone();v[:,pa]=c*up-s*uc;v[:,ch]=s*up+c*uc
  J[:,j,pa]=-alpha*v[:,ch];J[:,j,ch]=alpha*v[:,pa]
 K=w@w.transpose(1,2);S=(u@u.transpose(1,2)).square();yn=(K@S*K).sum(2).reshape(-1)
 return (8*J.square().sum(2)*yn[:,None]).clamp(min=1e-16)
best=1e9;began=time.time();accepted=0;stale=0
for step in range(a.steps):
 with torch.no_grad():w=project(decode(),w,a.cycles if step else max(100,a.cycles))
 theta.grad=None;u=decode();e2=loss(u,w);e2.backward();g=theta.grad.detach();value=float(e2.detach());error=np.sqrt(max(0,value))
 with torch.no_grad():
  if error<best:
   best=error;np.savez(a.out,U=u.cpu().numpy(),w=w.cpu().numpy(),W=np.zeros((R,C)),angle_words=theta.cpu().numpy().astype(np.int64)%(1<<a.brot),tree_edges=tree,brot=a.brot,err2=error)
  h=curvature(u,w);order=torch.argsort((-g.abs()+.5*h).flatten())[:max(512,16*a.moves)]
  choices=[];seen=set()
  for ix in order.cpu().tolist():
   row=ix//T
   if row not in seen:choices.append(ix);seen.add(row)
   if len(choices)>=a.moves:break
  groups=[choices[:k] for k in sorted(set([a.moves,8,4,2,1]),reverse=True) if k<=len(choices)]
  groups.extend([[int(ix)] for ix in order[1:65].cpu().tolist()])
  changed=False
  for ids in groups:
   ix=torch.tensor(ids,device='cuda');direction=-g.flatten()[ix].sign();flat=theta.view(-1);flat[ix]+=direction
   nv=float(loss(decode(),w))
   if nv<value-1e-12:changed=True;accepted+=len(ids);break
   flat[ix]-=direction
  stale=0 if changed else stale+1
 if step%5==0 or stale:print(json.dumps(dict(step=step,error_2=error,best=best,accepted=accepted,changed=changed,seconds=time.time()-began)),flush=True)
 if stale>=3:break
