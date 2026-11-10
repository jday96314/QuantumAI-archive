"""Integer Givens descent alternating with PSD block variable projection."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
from lattice_qubo import solve as solve_qubo

p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--cycles',type=int,default=12)
p.add_argument('--moves',type=int,default=16);p.add_argument('--brot',type=int,default=8)
p.add_argument('--candidates',type=int,default=384);p.add_argument('--trial-cycles',type=int,default=4)
a=p.parse_args();torch.set_num_threads(2);torch.backends.cuda.matmul.allow_tf32=False
root=(Path(__file__).resolve().parent/'../../3').resolve();z=dict(np.load(a.init));R,B,N=z['U'].shape;C=z['w'].shape[2]
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
 cross=2*(J.reshape(R,B*T,N)@u.transpose(1,2))*(u@u.transpose(1,2))[:,:,None,:].expand(R,B,T,B).reshape(R,B*T,B)
 L=torch.linalg.cholesky(S+1e-11*torch.eye(B,**d))
 solved=torch.cholesky_solve(cross.transpose(1,2),L).transpose(1,2)
 projected=(cross*solved).sum(2).reshape(R*B,T)
 raw=2*J.square().sum(2);h=(4*((raw-projected).clamp(min=0)+.03*raw)*yn[:,None]).clamp(min=1e-16)
 return h,J,L

@torch.no_grad()
def proposals(u,w,g,h,J,L):
 barrier=.5*h-g.abs();flat=barrier.flatten()
 ids=torch.unique(torch.cat([flat.argsort()[:a.candidates],(barrier/h).flatten().argsort()[:a.candidates//2],(g.square()/h).flatten().argsort(descending=True)[:a.candidates//2]]))
 rows=ids//T;angles=ids%T;ensembles=rows//B;localrows=rows%B;M=len(ids)
 deriv=J.reshape(R*B,T,N)[rows,angles];uv=u.reshape(R*B,N)[rows]
 # Residualize each feature tangent against its ensemble's dictionary.
 D=(deriv[:,sp['ii']]*uv[:,sp['jj']]+uv[:,sp['ii']]*deriv[:,sp['jj']])*fac
 P=features(u);X=w.transpose(1,2)@P
 Y=(w@X)[ensembles,localrows]
 for r in range(R):
  take=torch.nonzero(ensembles==r,as_tuple=True)[0]
  cr=D[take]@P[r].T
  D[take]-=torch.cholesky_solve(cr.T,L[r]).T@P[r]
 DY=D@Y.T;H=4*((D@D.T)*(Y@Y.T)+DY*DY.T)
 signs=-g.flatten()[ids].sign();base=(H*signs[:,None]*signs[None,:]).cpu().numpy()
 diag=h.flatten()[ids].cpu().numpy();np.fill_diagonal(base,diag)
 model=np.block([[base,-base],[-base,base]])
 gg=g.flatten()[ids].abs().cpu().numpy();unary=np.r_[.5*diag-gg,.5*diag+gg]
 selected,energy=solve_qubo(model,unary,cap=2*a.moves,seed=91335+step)
 spin=np.zeros(M,dtype=int)
 for ix in selected:spin[ix%M]+=1 if ix<M else -1
 chosen=np.flatnonzero(spin);rawids=ids.cpu().numpy();rawsigns=signs.cpu().numpy()
 pairs=[(int(rawids[i]),float(spin[i]*rawsigns[i])) for i in chosen]
 order=flat.argsort()[:max(512,8*a.moves)]
 greedy=[(int(ix),float(-g.flatten()[ix].sign())) for ix in order[:a.moves]]
 return pairs,greedy,energy
best=1e9;began=time.time();accepted=0;stale=0
for step in range(a.steps):
 with torch.no_grad():w=project(decode(),w,a.cycles if step else max(100,a.cycles))
 theta.grad=None;u=decode();e2=loss(u,w);e2.backward();g=theta.grad.detach();value=float(e2.detach());error=np.sqrt(max(0,value))
 with torch.no_grad():
  if error<best:
   best=error;np.savez(a.out,U=u.cpu().numpy(),w=w.cpu().numpy(),W=np.zeros((R,C)),angle_words=theta.cpu().numpy().astype(np.int64)%(1<<a.brot),tree_edges=tree,brot=a.brot,err2=error)
  h,J,L=curvature(u,w);pairs,greedy,energy=proposals(u,w,g,h,J,L)
  sizes=sorted(set([len(pairs),a.moves,128,64,32,16,8,4,2,1]),reverse=True)
  groups=[pairs[:k] for k in sizes if 0<k<=len(pairs)]
  groups += [greedy[:k] for k in sizes if 0<k<=len(greedy)]
  changed=False
  for group in groups:
   ix=torch.tensor([x[0] for x in group],device='cuda');direction=torch.tensor([x[1] for x in group],**d);flat=theta.view(-1);flat[ix]+=direction
   newu=decode();neww=project(newu,w,a.trial_cycles);nv=float(loss(newu,neww))
   if nv<value-1e-12:changed=True;accepted+=len(group);w=neww;break
   flat[ix]-=direction
  stale=0 if changed else stale+1
 if step%5==0 or stale:print(json.dumps(dict(step=step,error_2=error,best=best,accepted=accepted,changed=changed,proposal_energy=energy,seconds=time.time()-began)),flush=True)
 if stale>=3:break
