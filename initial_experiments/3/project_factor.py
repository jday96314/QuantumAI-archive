"""Project a fixed dictionary onto a requested physical factor rank."""
import argparse
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True);p.add_argument('--columns',type=int,required=True)
a=p.parse_args();torch.set_num_threads(1);z=dict(np.load(a.fit));root=Path(__file__).resolve().parent;sp=np.load(root/'parameters/reference_spectrum.npz')
u=z['U'][0];B=len(u);d=dict(device='cuda',dtype=torch.float64);P=torch.tensor(u[:,sp['ii']]*u[:,sp['jj']]*sp['fac'],**d)
G=torch.tensor((sp['eigenvectors']*np.maximum(sp['eigenvalues'],0))@sp['eigenvectors'].T,**d);S=P@P.T;H=P@G@P.T;L=torch.linalg.cholesky(S+1e-5*torch.eye(B,**d))
hh=torch.linalg.solve_triangular(L,H,upper=False);hh=torch.linalg.solve_triangular(L,hh.T,upper=False).T;e,v=torch.linalg.eigh(.5*(hh+hh.T))
w=torch.linalg.solve_triangular(L.T,v[:,-a.columns:]*e[-a.columns:].clamp(min=0).sqrt(),upper=True);K=w@w.T;KS=K@S
error=float((G.square().sum()-2*(K*H).sum()+(KS*KS.T).sum()).clamp(min=0).sqrt())
for key in ['K','onebody_angle_words','onebody_brot','error_1_optimized']:z.pop(key,None)
z.update(w=w.cpu().numpy()[None,:,:],W=np.zeros((1,a.columns)),err2=error);np.savez(a.out,**z);print('error_2',error,flush=True)
