"""Measure tensor error versus per-ensemble column rank, with fixed vectors."""
import numpy as np,time,json
from pathlib import Path
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
HERE=Path(__file__).resolve().parent
z=np.load(HERE/'candidates/chem35c50_rot/submission.npz')
init=np.load(HERE/'parameters/chem35_init.npz')
U=reconstruct_u_vectors(z['angle_words'][76:]/2048,z['tree_edges']).reshape(35,60,76)
ii,jj=np.tril_indices(76);fac=np.where(ii==jj,1,np.sqrt(2.))
P=(U[:,:,ii]*U[:,:,jj])*fac
X=np.einsum('rbc,rbp->rcp',init['w'],P)
val,vec=np.linalg.eigh(X@X.transpose(0,2,1));vec=vec[:,:,::-1]
xp=np.einsum('rcd,rcp->rdp',vec,X)
ref=np.load(HERE/'parameters/chem_reference.npz')['G']
ix=ii*76+jj;G=ref[ix[:,None],ix[None,:]].astype(float)*fac[:,None]*fac[None,:]
for C in [60,50,45,40,35,30,25,20,15,10,6,4]:
    y=xp[:,:C].reshape(-1,len(ii));d=y.T@y-G
    print(json.dumps(dict(C=C,error_2=float(np.linalg.norm(d)))),flush=True)
