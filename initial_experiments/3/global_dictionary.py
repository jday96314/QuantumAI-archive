"""Greedy ISDF-style rank-one dictionary from all current vector rows."""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--basis',type=int,default=500);p.add_argument('--out',required=True);p.add_argument('--fit')
a=p.parse_args();HERE=Path(__file__).resolve().parent
if a.fit:U=np.load(a.fit)['U'].reshape(-1,76)
else:
    z=np.load(HERE/'candidates/chem35c45_rot.npz');U=reconstruct_u_vectors(z['angle_words'][76:]/2048,z['tree_edges'])
sp=np.load(HERE/'parameters/reference_spectrum.npz');ev=sp['eigenvalues'];V=sp['eigenvectors'];ii,jj=sp['ii'],sp['jj'];fac=sp['fac']
P=U[:,ii]*U[:,jj]*fac;F=(V[:,:600]*np.sqrt(np.maximum(ev[:600],0))).T
gram=(U@U.T)**2;cross=P@F.T;g=gram.copy();cr=cross.copy();chosen=[];start=time.time()
for k in range(a.basis):
    diag=np.maximum(g.diagonal(),1e-20);score=(cr*cr).sum(1)/diag
    score[chosen]=-1;j=int(score.argmax());chosen.append(j)
    v=g[:,j]/np.sqrt(diag[j]);h=cr[j]/np.sqrt(diag[j]);g-=v[:,None]*v[None,:];cr-=v[:,None]*h[None,:]
    if k%100==0:print('select',k,time.time()-start,flush=True)
u=U[chosen];P=P[chosen];s=P@P.T
se,sv=eigh(s);se=np.maximum(se,1e-10);iroot=(sv/np.sqrt(se))@sv.T
Q=iroot@P;G=(V*ev)@V.T;H=Q@G@Q.T
e,v=eigh(H);e=np.maximum(e,0);order=np.argsort(e)[::-1];e=e[order];v=v[:,order]
w=(iroot@v)*np.sqrt(e);K=w@w.T
print(json.dumps(dict(B=a.basis,minimum_error_2=float(np.sqrt(max(0,np.square(ev).sum()-np.square(e).sum()))),
                      lambda_gram_bound=float(.25*abs(K).sum()),lambda_eigen_factor=float(.25*(abs(w).sum(0)**2).sum()))),flush=True)
np.savez(a.out,U=u[None,:,:],w=w[None,:,:],W=np.zeros((1,a.basis)),chosen=np.array(chosen))
