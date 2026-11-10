"""Exact rank lower bounds for the fixed shifted two-body tensor."""
from pathlib import Path
import numpy as np,json,time
from scipy.linalg import eigh
HERE=Path(__file__).resolve().parent
g=np.load(HERE/'parameters/chem_reference.npz')['G'];ii,jj=np.tril_indices(76)
ix=ii*76+jj;fac=np.where(ii==jj,1,np.sqrt(2.))
g=g[ix[:,None],ix[None,:]].astype(float)*fac[:,None]*fac[None,:]
t=time.time();ev,v=eigh(g,driver='evd');order=np.argsort(ev)[::-1];ev=ev[order];v=v[:,order]
np.savez(HERE/'parameters/reference_spectrum.npz',eigenvalues=ev,eigenvectors=v,ii=ii,jj=jj,fac=fac)
for rank in [100,150,200,250,292,300,350,400,500,600,700,900,1200,1500]:
    print(json.dumps(dict(rank=rank,min_error_2=float(np.linalg.norm(ev[rank:])))),flush=True)
print('seconds',time.time()-t,flush=True)
