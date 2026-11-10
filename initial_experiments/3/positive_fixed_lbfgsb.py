"""Bound-constrained factor fitting for a fixed quantized dictionary."""
import argparse, io, json, time, zipfile
from pathlib import Path
import numpy as np
import torch
from scipy.optimize import minimize

p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=2000);a=p.parse_args();torch.set_num_threads(1)
root=Path(__file__).resolve().parent
for _ in range(80):
 try:
  with np.load(io.BytesIO(Path(a.fit).read_bytes())) as z:z={k:z[k] for k in z.files}
  break
 except (OSError,ValueError,zipfile.BadZipFile,EOFError):time.sleep(.05)
sp=np.load(root/'parameters/reference_spectrum.npz');u=z['U'][0];w0=z['w'][0];B,C=w0.shape
P=u[:,sp['ii']]*u[:,sp['jj']]*sp['fac'];G=(sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T
d=dict(device='cuda',dtype=torch.float64);pt=torch.tensor(P,**d);gt=torch.tensor(G,**d)
S=pt@pt.T;H=pt@gt@pt.T;g2=gt.square().sum()
identity=(sp['ii']==sp['jj']).astype(float);t=sp['eigenvectors'].T@identity
trace_norm=float(np.maximum(sp['eigenvalues'],0)@(t*t));best=float('inf');began=time.time();calls=0
def fun(x):
 global best,calls
 w=torch.tensor(x.reshape(B,C),requires_grad=True,**d)
 ww=w*torch.sqrt(trace_norm/w.sum(0).square().sum())
 Z=S@ww;small=ww.T@Z
 e2=g2-2*(ww*(H@ww)).sum()+small.square().sum()
 e2.backward();err=float(e2.detach().clamp(min=0).sqrt());calls+=1
 if err<best:
  best=err;out=dict(z);out.update(w=ww.detach().cpu().numpy()[None,:,:],W=np.zeros((1,C)),err2=err)
  np.savez(a.out,**out)
 if calls%20==0 or calls==1:print(json.dumps(dict(calls=calls,B=B,error_2=err,best=best,seconds=time.time()-began)),flush=True)
 return float(e2.detach()),w.grad.detach().cpu().numpy().ravel()
r=minimize(fun,w0.ravel(),jac=True,bounds=[(0,None)]*(B*C),method='L-BFGS-B',options=dict(maxiter=a.steps,maxcor=15,ftol=1e-15,gtol=1e-10,maxls=30))
print(str(r.message),flush=True)
