"""Monotone integer-angle moves alternating with nonnegative weight fitting."""
import argparse,io,json,time,zipfile
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=4000);p.add_argument('--weight-steps',type=int,default=10);p.add_argument('--weight-lr',type=float,default=1e-5);p.add_argument('--moves',type=int,default=32)
p.add_argument('--projected',action='store_true');p.add_argument('--trial-weight-steps',type=int,default=0)
a=p.parse_args();torch.set_num_threads(1);root=Path(__file__).resolve().parent
for _ in range(60):
 try:
  with np.load(io.BytesIO(Path(a.fit).read_bytes())) as z:z={k:z[k] for k in z.files}
  break
 except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
R,B,N=z['U'].shape;assert R==1;C=z['w'].shape[2];bits=int(z.get('brot',8));tree=z['tree_edges'];sp=np.load(root/'parameters/reference_spectrum.npz')
d=dict(device='cuda',dtype=torch.float32);G=torch.tensor((sp['eigenvectors']*np.maximum(sp['eigenvalues'],0))@sp['eigenvectors'].T,**d)
ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
identity=(sp['ii']==sp['jj']).astype(float);tr=sp['eigenvectors'].T@identity;trace_norm=float(np.maximum(sp['eigenvalues'],0)@(tr*tr))
theta=torch.tensor(z.get('angle_words',encode(z['U'][0],tree,bits)).reshape(B,75),**d)
w=torch.tensor(np.maximum(z['w'][0],0),requires_grad=True,**d);opt=torch.optim.Adam([w],lr=a.weight_lr)
start=torch.zeros((B,N),**d);start[:,0]=1;alpha=np.pi/(1<<(bits-1));stages=[]
for j,(pa,ch) in enumerate(tree):
 pa,ch=int(pa),int(ch)
 if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
 stages[-1].append((j,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
def decode():
 u=start
 for ts,ps,cs in stages:
  ang=theta[:,ts]*alpha;c,s=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs];new=u.clone();new[:,ps]=c*up-s*uc;new[:,cs]=s*up+c*uc;u=new
 return u
def weights():return w*torch.sqrt(trace_norm/w.sum(0).square().sum().clamp(min=1e-15))
def loss(P):
 X=weights().T@P;return (X.T@X-G).square().sum()
def features(u):return u[:,ii]*u[:,jj]*fac
@torch.no_grad()
def curvature(u,g):
 v=start.clone();J=torch.zeros((B,75,N),**d)
 for j,(pa,ch) in enumerate(tree):
  pa,ch=int(pa),int(ch);ang=theta[:,j]*alpha;c,s=ang.cos(),ang.sin();jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone()
  J[:,:,pa]=c[:,None]*jp-s[:,None]*jc;J[:,:,ch]=s[:,None]*jp+c[:,None]*jc
  up,uc=v[:,pa].clone(),v[:,ch].clone();v[:,pa]=c*up-s*uc;v[:,ch]=s*up+c*uc;J[:,j,pa]=-alpha*v[:,ch];J[:,j,ch]=alpha*v[:,pa]
 ww=weights();K=ww@ww.T;S=(u@u.T).square();yn=(K@S*K).sum(1)
 h=(8*J.square().sum(2)*yn[:,None]).clamp(min=1e-16)
 if a.projected:
  barrier=.5*h-g.abs();energy=g.square()/h;ids=torch.unique(torch.cat([barrier.flatten().argsort()[:512],energy.flatten().argsort(descending=True)[:512]]))
  rows=ids//75;j=J.reshape(B*75,N)[ids].double();ud=u.double();cross=2*(j@ud.T)*(ud[rows]@ud.T)
  L=torch.linalg.cholesky((ud@ud.T).square()+1e-5*torch.eye(B,device='cuda',dtype=torch.float64))
  projected=(cross*torch.cholesky_solve(cross.T,L).T).sum(1);raw=2*j.square().sum(1)
  h.view(-1)[ids]=(4*((raw-projected).clamp(min=0)+.1*raw)*yn[rows].double()).float().clamp(min=1e-16)
 return h
def refit(P,count):
 with torch.no_grad():bestloss=float(loss(P));bestw=w.detach().clone()
 for _ in range(count):
  with torch.enable_grad():opt.zero_grad();error=loss(P);error.backward();opt.step()
  with torch.no_grad():
   w.clamp_(min=0);value=float(loss(P))
   if value<bestloss:bestloss=value;bestw=w.detach().clone()
 with torch.no_grad():w.copy_(bestw)
 return bestloss
best=1e9;began=time.time();accepted=0;stale=0
for step in range(a.steps):
 theta.requires_grad_(False)
 with torch.no_grad():u=decode();P=features(u);base=float(loss(P));bestw=w.detach().clone();bestweight=base
 bestweight=refit(P,a.weight_steps)
 theta.requires_grad_(True);u=decode();P=features(u);e2=loss(P);g=torch.autograd.grad(e2,theta)[0].detach();value=float(e2.detach());err=np.sqrt(value)
 with torch.no_grad():
  if err<best:
   best=err;np.savez(a.out,U=u.cpu().numpy().astype(float)[None,:,:],w=weights().cpu().numpy().astype(float)[None,:,:],W=np.zeros((1,C)),tree_edges=tree,brot=bits,angle_words=theta.cpu().numpy().astype(np.int64)%(1<<bits),err2=err)
  h=curvature(u,g);scores=.5*h-g.abs();order=scores.flatten().argsort()[:512];choices=[];seen=set()
  for ix in order.cpu().tolist():
   row=ix//75
   if row not in seen:choices.append(ix);seen.add(row)
   if len(choices)>=a.moves:break
  sizes=sorted(set([a.moves>>k for k in range(a.moves.bit_length())]+[16,8,4,2,1]),reverse=True)
  groups=[choices[:k] for k in sizes if k<=len(choices)]+[[int(ix)] for ix in order[1:(9 if a.trial_weight_steps else 65)].cpu().tolist()]
  changed=False
  for ids in groups:
   ix=torch.tensor(ids,device='cuda');direction=-g.flatten()[ix].sign();flat=theta.view(-1);flat[ix]+=direction
   newP=features(decode());trial=float(loss(newP));savedw=None;savedstate=None
   if trial>=value-1e-11 and a.trial_weight_steps:
    savedw=w.detach().clone();savedstate={k:v.clone() if torch.is_tensor(v) else v for k,v in opt.state[w].items()}
    trial=refit(newP,a.trial_weight_steps)
   if trial<value-1e-11:changed=True;accepted+=len(ids);break
   if savedw is not None:w.copy_(savedw);opt.state[w].update(savedstate)
   flat[ix]-=direction
  stale=0 if changed or base-bestweight>1e-10 else stale+1
 if step%10==0 or stale:print(json.dumps(dict(step=step,error_2=err,best=best,accepted=accepted,changed=changed,seconds=time.time()-began)),flush=True)
 if stale>=10:break
