"""Joint lattice continuation and constrained normalization refinement."""
import argparse,io,json,time,zipfile
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=24000);p.add_argument('--angle-lr',type=float,default=.002);p.add_argument('--weight-lr',type=float,default=2e-6)
p.add_argument('--cap',type=float,default=.01998);p.add_argument('--strength',type=float,default=100.)
a=p.parse_args();torch.set_num_threads(1);root=Path(__file__).resolve().parent
for _ in range(60):
 try:
  with np.load(io.BytesIO(Path(a.fit).read_bytes())) as z:z={k:z[k] for k in z.files}
  break
 except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
R,B,N=z['U'].shape;assert R==1;C=z['w'].shape[2];bits=int(z.get('brot',8));tree=z['tree_edges'];alpha=np.pi/(1<<(bits-1))
d=dict(device='cuda',dtype=torch.float64);theta=torch.tensor(z.get('angle_words',encode(z['U'][0],tree,bits)).reshape(B,75),requires_grad=True,**d)
w=torch.tensor(z['w'][0],requires_grad=True,**d);off=torch.tensor(z['W'][0],requires_grad=True,**d)
sp=np.load(root/'parameters/reference_spectrum.npz');ref=np.load(root/'parameters/chem_reference.npz');ev,V=sp['eigenvalues'],sp['eigenvectors']
G=torch.tensor((V*ev)@V.T,device='cuda',dtype=torch.float32);ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
h0=torch.tensor(ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(N),**d)
constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1']);start=torch.zeros((B,N),**d);start[:,0]=1
stages=[]
for j,(pa,ch) in enumerate(tree):
 pa,ch=int(pa),int(ch)
 if not stages or any(pa in (x,y) or ch in (x,y) for _,x,y in stages[-1]):stages.append([])
 stages[-1].append((j,pa,ch))
stages=[tuple(torch.tensor(x,device='cuda') for x in zip(*s)) for s in stages]
def decode(t):
 u=start
 for ts,ps,cs in stages:
  ang=t[:,ts]*alpha;c,s=ang.cos(),ang.sin();up,uc=u[:,ps],u[:,cs];new=u.clone();new[:,ps]=c*up-s*uc;new[:,cs]=s*up+c*uc;u=new
 return u
def vals(t):
 u=decode(t);P=u[:,ii]*u[:,jj]*fac;X=w.T@P;xf=X.float();e2=(xf.T@xf-G).square().sum().double()
 packed=(off@X)/fac;hc=torch.zeros((N,N),**d);hc[ii,jj]=packed;hc=hc+hc.T-torch.diag(hc.diag());h=h0-hc
 e=torch.linalg.eigvalsh(.5*(h+h.T));lh1=torch.sqrt(e.square()+1e-14).sum();gap=constant+lh1-e.sum()+.5*off.square().sum()
 lam=lh1+.25*(torch.sqrt(w.square()+1e-14).sum(0)+torch.sqrt((off+w.sum(0)).square()+1e-14)).square().sum()
 return e2,gap,lam,gap*(2*lam-gap),u
with torch.no_grad():
 u=decode(theta);K=w@w.T;S=(u@u.T).square();yn=(K@S*K).sum(1);v=start.clone();J=torch.zeros((B,75,N),**d)
 for j,(pa,ch) in enumerate(tree):
  pa,ch=int(pa),int(ch);ang=theta[:,j]*alpha;c,s=ang.cos(),ang.sin();jp,jc=J[:,:,pa].clone(),J[:,:,ch].clone()
  J[:,:,pa]=c[:,None]*jp-s[:,None]*jc;J[:,:,ch]=s[:,None]*jp+c[:,None]*jc
  up,uc=v[:,pa].clone(),v[:,ch].clone();v[:,pa]=c*up-s*uc;v[:,ch]=s*up+c*uc;J[:,j,pa]=-alpha*v[:,ch];J[:,j,ch]=alpha*v[:,pa]
 curvature=(8*J.square().sum(2)*yn[:,None]).clamp(min=1e-16)
opt=torch.optim.Adam([{'params':[theta],'lr':a.angle_lr},{'params':[w],'lr':a.weight_lr},{'params':[off],'lr':3*a.weight_lr}])
mue=float(z.get('lagrange_error',.5));mug=float(z.get('lagrange_gap',.6));rho=2.;best=float('inf');began=time.time()
def checkpoint(step):
 global best
 with torch.no_grad():
  words=theta.round();ee,gg,ll,ss,uu=vals(words);error=float(ee.sqrt());gap=float(gg);value=float(ss.clamp(min=0).sqrt())
  if error<=.019995 and gap>=3.9151 and value<best:
   best=value;out=dict(U=uu.cpu().numpy()[None,:,:],w=w.cpu().numpy()[None,:,:],W=off.cpu().numpy()[None,:],tree_edges=tree,brot=bits,angle_words=words.cpu().numpy().astype(np.int64)%(1<<bits),err2=error,lagrange_error=mue,lagrange_gap=mug)
   np.savez(a.out,**out)
 return error,gap,value
checkpoint(-1)
for step in range(a.steps):
 opt.zero_grad();ee,gg,ll,ss,u=vals(theta);ge=ee/(a.cap*a.cap)-1;gp=(3.9155-gg)/3.9155
 def penalty(g,mu):return (torch.relu(mu+rho*g).square()-mu*mu)/(2*rho)
 strength=10**(-3+(np.log10(a.strength)+3)*step/max(1,a.steps-1));lattice=(curvature*(1-torch.cos(2*np.pi*theta))/(4*np.pi*np.pi)).sum()/(a.cap*a.cap)
 objective=ss/3000+penalty(ge,mue)+penalty(gp,mug)+strength*lattice;objective.backward();opt.step()
 decay=.1+.9*.5*(1+np.cos(np.pi*step/a.steps));opt.param_groups[0]['lr']=a.angle_lr*decay;opt.param_groups[1]['lr']=a.weight_lr*decay;opt.param_groups[2]['lr']=3*a.weight_lr*decay
 if step%20==0:mue=max(0,mue+rho*float(ge.detach()));mug=max(0,mug+rho*float(gp.detach()))
 if step%100==0 or step==a.steps-1:
  error,gap,value=checkpoint(step)
  with torch.no_grad():
   np.savez(str(Path(a.out).with_suffix(''))+'_soft.npz',U=decode(theta).cpu().numpy()[None,:,:],w=w.cpu().numpy()[None,:,:],W=off.cpu().numpy()[None,:],tree_edges=tree,err2=float(ee.detach().sqrt()))
  print(json.dumps(dict(step=step,soft_error=float(ee.detach().sqrt()),hard_error=error,hard_gap=gap,hard_lambda=value,best=best,strength=strength,seconds=time.time()-began)),flush=True)
