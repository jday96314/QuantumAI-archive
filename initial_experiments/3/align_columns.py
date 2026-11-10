"""Preserve a fit while rotating its factor columns toward a sparse reference."""
import argparse,json
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--reference',required=True);p.add_argument('--out',required=True)
p.add_argument('--weights-only',action='store_true')
p.add_argument('--physical',action='store_true',help='align represented matrices across different dictionaries')
a=p.parse_args();z=dict(np.load(a.fit));ref=np.load(a.reference);R,B,C=z['w'].shape
if not a.physical:assert ref['w'].shape==z['w'].shape
else:assert ref['w'].shape[0]==R and ref['w'].shape[2]>=C
for r in range(R):
 A=np.concatenate([z['w'][r],(z['W'][r]+z['w'][r].sum(0))[None,:]],0)
 if a.physical:
  ur=ref['U'][r];un=z['U'][r];wr=ref['w'][r];ii,jj=np.tril_indices(un.shape[1]);fac=np.where(ii==jj,1,np.sqrt(2))
  xn=z['w'][r].T@(un[:,ii]*un[:,jj]*fac);xr=wr.T@(ur[:,ii]*ur[:,jj]*fac)
  chosen=np.argsort(np.square(xr).sum(1))[::-1][:C];cross=xn@xr[chosen].T
 else:
  target=np.concatenate([ref['w'][r],(ref['W'][r]+ref['w'][r].sum(0))[None,:]],0)
  cross=(A[:B].T@target[:B]) if a.weights_only else (A.T@target)
 left,_,right=np.linalg.svd(cross,full_matrices=False);q=left@right;rot=A@q
 assert np.allclose(A@A.T,rot@rot.T,rtol=1e-10,atol=1e-10)
 z['w'][r]=rot[:B];z['W'][r]=rot[B]-rot[:B].sum(0)
 print(json.dumps(dict(r=r,before=float(.25*np.square(abs(A).sum(0)).sum()),after=float(.25*np.square(abs(rot).sum(0)).sum()))),flush=True)
np.savez(a.out,**z)
