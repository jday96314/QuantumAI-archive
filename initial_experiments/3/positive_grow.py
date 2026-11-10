"""Grow a nonnegative quantized dictionary using the physical loss gradient."""
import argparse, io, json, time, zipfile
from pathlib import Path
import numpy as np
import torch
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--reference',action='append',required=True)
p.add_argument('--out',required=True);p.add_argument('--basis',type=int,default=1000);p.add_argument('--refit-steps',type=int,default=30)
p.add_argument('--lr',type=float,default=1e-5);a=p.parse_args();torch.set_num_threads(1)
root=Path(__file__).resolve().parent
for _ in range(80):
 try:
  with np.load(io.BytesIO(Path(a.fit).read_bytes())) as z:z={k:z[k] for k in z.files}
  break
 except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
bits=int(z['brot']);tree=z['tree_edges'];B=z['U'].shape[1];C=z['w'].shape[2]
rng=np.random.default_rng(39071);tables=[]
for ref in a.reference:
 u=np.load(ref)['U'][0].copy();angles=np.zeros((len(u),75))
 for j,(pa,ch) in reversed(list(enumerate(tree))):
  angles[:,j]=np.arctan2(u[:,ch],u[:,pa]);u[:,pa]=np.hypot(u[:,pa],u[:,ch]);u[:,ch]=0
 units=angles/np.pi*(1<<(bits-1));lo=np.floor(units);frac=units-lo
 tables.append(np.rint(units))
 tables.extend(lo+(frac>t) for t in [.2,.4,.6,.8])
 tables.extend(lo+(frac>rng.random(frac.shape)) for _ in range(2))
words=np.unique(np.concatenate(tables).astype(np.int64)%(1<<bits),axis=0)
oldwords=z['angle_words'].reshape(B,75);oldset={tuple(row) for row in oldwords}
words=np.array([row for row in words if tuple(row) not in oldset])
uc=reconstruct_u_vectors(words/(1<<(bits-1)),tree)
sp=np.load(root/'parameters/reference_spectrum.npz');ii,jj,fac=sp['ii'],sp['jj'],sp['fac']
d=dict(device='cuda',dtype=torch.float32)
G=torch.tensor((sp['eigenvectors']*np.maximum(sp['eigenvalues'],0))@sp['eigenvectors'].T,**d)
P=torch.tensor(z['U'][0][:,ii]*z['U'][0][:,jj]*fac,**d);Pc=torch.tensor(uc[:,ii]*uc[:,jj]*fac,**d)
identity=(ii==jj).astype(float);t=sp['eigenvectors'].T@identity;trace_norm=float(np.maximum(sp['eigenvalues'],0)@(t*t))
w=torch.tensor(z['w'][0],**d);chosen=[];began=time.time();best=float('inf')
def normalized():return w*torch.sqrt(trace_norm/w.sum(0).square().sum())
def loss():
 X=normalized().T@P
 return (X.T@X-G).square().sum()
def save(err):
 out=dict(U=np.r_[z['U'][0],uc[chosen]][None,:,:],w=normalized().detach().cpu().numpy().astype(float)[None,:,:],W=np.zeros((1,C)),tree_edges=tree,brot=bits,angle_words=np.r_[oldwords,words[chosen]],err2=err)
 np.savez(a.out,**out)
for growth in range(a.basis-B+1):
 with torch.no_grad():
  ww=normalized();X=ww.T@P;pred=X.T@X;res=pred-G;err=float(res.square().sum().sqrt())
  save(err)
  print(json.dumps(dict(B=len(P),error_2=err,lambda_h2_no_offset=float(.25*ww.sum(0).square().sum()),seconds=time.time()-began)),flush=True)
  if growth==a.basis-B or err<.01985:break
  grad=4*(Pc@(res@X.T)) - 4*(res*pred).sum()/trace_norm*ww.sum(0)[None,:]
  curv=4*(X.square().sum(1)[None,:]+(Pc@X.T).square()).clamp(min=1e-8)
  scores=(grad.clamp(max=0).square()/curv).sum(1)
  if chosen:scores[torch.tensor(chosen,device='cuda')]=-1
  ix=int(scores.argmax());chosen.append(ix)
  initial=(-grad[ix]).clamp(min=0);initial*=1e-5/initial.max().clamp(min=1e-15)
  w=torch.cat([ww,initial[None,:]],0).detach().requires_grad_(True);P=torch.cat([P,Pc[ix:ix+1]],0)
 opt=torch.optim.Adam([w],lr=a.lr)
 with torch.no_grad():bestloss=float(loss());bestw=w.detach().clone()
 for _ in range(a.refit_steps):
  opt.zero_grad();e=loss();e.backward();opt.step()
  with torch.no_grad():
   w.clamp_(min=0);value=float(loss())
   if value<bestloss:bestloss=value;bestw=w.detach().clone()
 with torch.no_grad():w.copy_(bestw)
