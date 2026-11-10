"""Add residual-directed global rank-one features to each THC ensemble."""
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True);p.add_argument('--basis',type=int,required=True)
a=p.parse_args();root=Path(__file__).resolve().parent;z=np.load(a.fit);u,w=z['U'],z['w'];R,B,C=w.shape
sp=np.load(root/'parameters/reference_spectrum.npz');ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];ev,V=sp['eigenvalues'],sp['eigenvectors'];G=(V*np.maximum(ev,0))@V.T
pool=np.r_[np.load(root/'parameters/global500_lbfgs_fast_init.npz')['U'][0],np.load(root/'parameters/positive725_valid_float.npz')['U'][0]]
Pc=pool[:,ii]*pool[:,jj]*fac;gram=Pc@Pc.T;GP=Pc@G
P=u[:,:,ii]*u[:,:,jj]*fac;X=np.einsum('rbc,rbp->rcp',w,P).reshape(R*C,-1);residual=GP-(Pc@X.T)@X
outs=[];ows=[]
for r in range(R):
    old=P[r];S=old@old.T;si=np.linalg.inv(S+1e-9*np.eye(B));cp=Pc@old.T
    work=gram-cp@si@cp.T;xr=X[r*C:(r+1)*C]
    cross=residual@xr.T-cp@si@(old@(G-X.T@X)@xr.T);chosen=[]
    for k in range(a.basis-B):
        diag=np.maximum(np.diag(work),1e-15);score=np.sum(cross*cross,axis=1)/diag;score[diag<1e-10]=-1
        if chosen:score[chosen]=-1
        j=int(score.argmax());chosen.append(j);v=work[:,j]/np.sqrt(diag[j]);h=cross[j]/np.sqrt(diag[j]);work-=np.outer(v,v);cross-=np.outer(v,h)
    outs.append(np.r_[u[r],pool[chosen]]);ows.append(np.r_[w[r],np.zeros((a.basis-B,C))]);print(json.dumps(dict(ensemble=r,added=len(chosen))),flush=True)
np.savez(a.out,U=np.array(outs),w=np.array(ows),W=np.zeros((R,C)))
