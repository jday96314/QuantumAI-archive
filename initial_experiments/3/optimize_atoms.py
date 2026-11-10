"""Optimize residual-directed atoms, round them, and verify single swaps."""
import argparse,io,json,time,zipfile
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--pool',required=True);p.add_argument('--out',required=True)
p.add_argument('--cycles',type=int,default=8);p.add_argument('--steps',type=int,default=400);p.add_argument('--candidates',type=int,default=128)
p.add_argument('--columns',type=int,default=340);p.add_argument('--lattice',action='store_true');p.add_argument('--angle-lr',type=float,default=.1)
a=p.parse_args();torch.set_num_threads(1);root=Path(__file__).resolve().parent
for _ in range(60):
 try:
  with np.load(io.BytesIO(Path(a.fit).read_bytes())) as z:z={k:z[k] for k in z.files}
  break
 except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
pool=np.load(a.pool);u=z['U'][0];words=z['angle_words'].reshape(-1,75);tree=z['tree_edges'];bits=int(z['brot']);assert np.array_equal(tree,pool['tree_edges'])
sp=np.load(root/'parameters/reference_spectrum.npz');ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];d=dict(device='cuda',dtype=torch.float64)
ev,V=sp['eigenvalues'],sp['eigenvectors'];ix=np.argsort(ev)[-a.columns:];Z=torch.tensor(V[:,ix]*np.sqrt(np.maximum(ev[ix],0)),**d)
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);const=G.square().sum();ti=torch.tensor(ii,device='cuda');tj=torch.tensor(jj,device='cuda');ff=torch.tensor(fac,**d)
tried=set();began=time.time();accepted=0;rng=np.random.default_rng(78102)
stages=[]
for j,(pa,ch) in enumerate(tree):
 pa,ch=int(pa),int(ch)
 if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
 stages[-1].append((j,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
def decode(theta):
 angles=(theta+(theta.round()-theta).detach())*(np.pi/(1<<(bits-1)))
 un=torch.zeros((len(theta),76),**d);un[:,0]=1
 for ts,ps,cs in stages:
  ang=angles[:,ts];c,s=ang.cos(),ang.sin();up,uc=un[:,ps],un[:,cs];new=un.clone()
  new[:,ps]=c*up-s*uc;new[:,cs]=s*up+c*uc;un=new
 return un
def feature(v):return v[:,ti]*v[:,tj]*ff
def rank_error(S,H):
 L=torch.linalg.cholesky(S+1e-5*torch.eye(len(S),**d));hh=torch.linalg.solve_triangular(L,H,upper=False);hh=torch.linalg.solve_triangular(L,hh.T,upper=False).T
 ee,vv=torch.linalg.eigh(.5*(hh+hh.T));f=torch.linalg.solve_triangular(L.T,vv[:,-a.columns:]*ee[-a.columns:].clamp(min=0).sqrt(),upper=True);K=f@f.T;KS=K@S
 return float((const-2*(K*H).sum()+(KS*KS.T).sum()).clamp(min=0).sqrt())
for cycle in range(a.cycles):
 P=feature(torch.tensor(u,**d));S=P@P.T;H=P@G@P.T;si=torch.linalg.inv(S+1e-10*torch.eye(len(S),**d));cap=si@H;K=cap@si
 loss=2*(cap*K.T).sum(1)/si.diag()-(K.diag()/si.diag()).square()
 for j,w in enumerate(words):
  if tuple(w) in tried:loss[j]=float('inf')
 remove=int(loss.argmin());tried.add(tuple(words[remove]));keep=torch.arange(len(u),device='cuda')!=remove
 P0=P[keep];S0=S[keep][:,keep];sib=(si-torch.outer(si[:,remove],si[remove,:])/si[remove,remove])[keep][:,keep]
 F0=P0@Z;H0=F0@F0.T;fullH0=H[keep][:,keep];original=rank_error(S,H)
 def gain(v,core=True):
  pp=feature(v);cross=P0@pp.T;t=sib@cross;diag=(pp.square().sum(1)-(cross*t).sum(0)).clamp(min=1e-8)
  if core:f=pp@Z;hc=F0@f.T;hd=f.square().sum(1);hh=H0
  else:pg=pp@G;hc=P0@pg.T;hd=(pp*pg).sum(1);hh=fullH0
  nr=hc-hh@t;hr=hd-2*(t*hc).sum(0)+(t*(hh@t)).sum(0)
  return (2*(nr*(sib@nr)).sum(0)/diag+(hr/diag).square()).clamp(min=1e-16)
 with torch.no_grad():
  initial=torch.tensor(pool['U'][0],**d);scores=gain(initial);order=scores.argsort(descending=True)[:a.candidates]
 v=initial[order].clone().requires_grad_(True);snapshots=[]
 if a.lattice:
  theta=torch.tensor(pool['angle_words'][order.cpu().numpy()],requires_grad=True,**d);opt=torch.optim.Adam([theta],lr=a.angle_lr)
 else:opt=torch.optim.Adam([v],lr=.002)
 bestscores=torch.full((len(v),),-1.,**d);bestvectors=v.detach().clone()
 for step in range(a.steps):
  opt.zero_grad();un=decode(theta) if a.lattice else v/v.norm(dim=1,keepdim=True);scores=gain(un);objective=-scores.log().mean();objective.backward();opt.step()
  with torch.no_grad():
   improved=scores>bestscores;bestvectors[improved]=un[improved];bestscores=torch.maximum(bestscores,scores)
  if step%100==0 or step==a.steps-1:
   snapshots.append(un.detach().cpu().numpy());print(json.dumps(dict(cycle=cycle,step=step,max_gain=float(scores.detach().max()),seconds=time.time()-began)),flush=True)
 snapshots.append(bestvectors.cpu().numpy());vectors=np.concatenate(snapshots);work=vectors.copy();theta=np.zeros((len(vectors),75))
 for j in reversed(range(75)):
  pa,ch=tree[j];theta[:,j]=np.arctan2(work[:,ch],work[:,pa]);work[:,pa]=np.hypot(work[:,pa],work[:,ch]);work[:,ch]=0
 units=theta/np.pi*(1<<(bits-1));lower=np.floor(units);fraction=units-lower
 tables=[np.rint(units)]+[lower+(fraction>t) for t in [.25,.5,.75]]+[lower+(fraction>rng.random(fraction.shape)) for _ in range(2)]
 cw=np.unique(np.r_[words[remove:remove+1],np.concatenate(tables).astype(np.int64)%(1<<bits)],axis=0)
 uc=reconstruct_u_vectors(cw/(1<<(bits-1)),tree)
 with torch.no_grad():
  pc=feature(torch.tensor(uc,**d));scores=gain(torch.tensor(uc,**d),False);order=scores.argsort(descending=True)[:12]
  best=original;bestj=None
  for j in order.cpu().tolist():
   cross=P0@pc[j];hc=P0@(G@pc[j]);ss=torch.empty_like(S);hh=torch.empty_like(H)
   ss[:-1,:-1]=S0;ss[:-1,-1]=cross;ss[-1,:-1]=cross;ss[-1,-1]=pc[j]@pc[j]
   hh[:-1,:-1]=fullH0;hh[:-1,-1]=hc;hh[-1,:-1]=hc;hh[-1,-1]=pc[j]@G@pc[j]
   error=rank_error(ss,hh)
   if error<best-1e-9:best=error;bestj=j
 if bestj is not None:
  u=np.r_[u[np.arange(len(u))!=remove],uc[bestj:bestj+1]];words=np.r_[words[np.arange(len(words))!=remove],cw[bestj:bestj+1]];accepted+=1
 np.savez(a.out,U=u[None,:,:],angle_words=words,tree_edges=tree,brot=bits,err2=best)
 print(json.dumps(dict(cycle=cycle,original_error=original,error_2=best,accepted=accepted,seconds=time.time()-began)),flush=True)
