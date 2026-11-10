"""Quantize a global dictionary and solve its PSD coefficient matrix exactly."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=12);p.add_argument('--columns',type=int,default=500)
p.add_argument('--ridge',type=float,default=0.)
a=p.parse_args();HERE=Path(__file__).resolve().parent
z=np.load(a.fit);U=z['U'][0];B=len(U);base=np.load(HERE/'candidates/chem35c45_rot.npz');tree=z['tree_edges'] if 'tree_edges' in z else base['tree_edges']
words=encode(U,tree,a.brot);U=reconstruct_u_vectors(words/(1<<(a.brot-1)),tree)
sp=np.load(HERE/'parameters/reference_spectrum.npz');ev,V=sp['eigenvalues'],sp['eigenvectors'];ii,jj=sp['ii'],sp['jj'];fac=sp['fac']
G=(V*ev)@V.T;P=U[:,ii]*U[:,jj]*fac;S=P@P.T
e,v=eigh(S);e=np.maximum(e,1e-12);ir=(v/np.sqrt(e))@v.T;Q=ir@P
H=Q@G@Q.T;e,v=eigh(H);order=np.argsort(e)[::-1];e=np.maximum(e[order],0);v=v[:,order]
if a.ridge:
    Sr=S+a.ridge*np.eye(B);Hraw=P@G@P.T
    Gpositive=(V*np.maximum(ev,0))@V.T
    Hpositive=P@Gpositive@P.T
    K=np.linalg.solve(Sr,np.linalg.solve(Sr,Hpositive).T).T;K=.5*(K+K.T)
    ek,vk=eigh(K);full=vk*np.sqrt(np.maximum(ek,0))
    X=full.T@P;ex,vx=eigh(X@X.T);ix=np.argsort(ex)[::-1]
    w_all=full@vx[:,ix]
    def ridge_error(C):
        ww=w_all[:,:C];kk=ww@ww.T;ks=kk@S
        return np.sqrt(max(0,np.square(ev).sum()-2*(kk*Hraw).sum()+(ks*ks.T).sum()))
for C in [300,350,400,450,500]:
    if C>B:continue
    error=(ridge_error(C) if a.ridge else np.sqrt(max(0,float(np.square(ev).sum()-np.square(e[:C]).sum()))))
    print(json.dumps(dict(brot=a.brot,columns=C,error_2=float(error))),flush=True)
C=min(a.columns,B);w=w_all[:,:C] if a.ridge else (ir@v[:,:C])*np.sqrt(e[:C]);K=w@w.T
np.savez(a.out,U=U[None,:,:],w=w[None,:,:],W=np.zeros((1,C)),angle_words=words,tree_edges=tree,brot=a.brot,
         err2=ridge_error(C) if a.ridge else np.sqrt(max(0,np.square(ev).sum()-np.square(e[:C]).sum())))
print('Gram LCU lower bound',float(.25*abs(K).sum()),flush=True)
