"""Backward feature selection by exact two-sided Frobenius projection loss."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import solve
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--basis',type=int,required=True);a=p.parse_args();HERE=Path(__file__).resolve().parent
z=np.load(a.fit);u=z['U'][0];sp=np.load(HERE/'parameters/reference_spectrum.npz')
ev,V=sp['eigenvalues'],sp['eigenvectors'];G=(V*np.maximum(ev,0))@V.T
P=u[:,sp['ii']]*u[:,sp['jj']]*sp['fac'];S=P@P.T;H=P@G@P.T
si=np.linalg.inv(S);chosen=np.arange(len(u));const=float(np.square(ev).sum());start=time.time()
while len(chosen)>a.basis:
 K=si@H@si;ksk=K@S@K;diag=np.diag(si)
 loss=2*np.diag(ksk)/diag-np.square(np.diag(K)/diag)
 ix=int(np.argmin(loss));keep=np.arange(len(chosen))!=ix
 si=(si-np.outer(si[:,ix],si[ix,:])/si[ix,ix])[np.ix_(keep,keep)]
 S=S[np.ix_(keep,keep)];H=H[np.ix_(keep,keep)];chosen=chosen[keep]
 if len(chosen)%10==0:
  cap=si@H;err=np.sqrt(max(0,const-(cap*cap.T).sum()))
  print(json.dumps(dict(B=len(chosen),error_2=float(err),last_loss=float(loss[ix]),seconds=time.time()-start)),flush=True)
out=dict(U=u[chosen][None,:,:],chosen=chosen)
if 'tree_edges' in z:out['tree_edges']=z['tree_edges']
np.savez(a.out,**out)
