"""Compress merged ensembles with a greedy rank-one dictionary selection."""
import argparse,json,time
from pathlib import Path
import numpy as np
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--group',type=int,default=3)
p.add_argument('--basis',type=int,default=128);p.add_argument('--columns',type=int,default=120);p.add_argument('--out',required=True)
a=p.parse_args();HERE=Path(__file__).resolve().parent;z=dict(np.load(a.source));R,B,C=z['W'].shape;n=76
U=reconstruct_u_vectors(z['angle_words'][n:]/2048,z['tree_edges']).reshape(R,B,n)
oldW=z['H2_I']-z['W'].sum(1)
ii,jj=np.tril_indices(n);fac=np.where(ii==jj,1,np.sqrt(2.));ix=ii*n+jj
Us=[];ws=[];Ws=[]
for start in range(0,R,a.group):
    count=min(a.group,R-start);u=U[start:start+count].reshape(-1,n)
    P=u[:,ii]*u[:,jj]*fac
    w=np.zeros((count*B,count*C))
    for k in range(count):w[k*B:(k+1)*B,k*C:(k+1)*C]=z['W'][start+k]
    of=oldW[start:start+count].reshape(-1)
    X=w.T@P
    if len(u)>a.basis:
        G=P@P.T;cross=P@X.T;work=G.copy();cr=cross.copy();piv=[]
        for step in range(a.basis):
            diag=np.maximum(np.diag(work),1e-30)
            scores=(cr*cr).sum(1)/diag
            scores[piv]=-1
            j=int(scores.argmax());piv.append(j)
            v=work[:,j]/np.sqrt(diag[j]);h=cr[j]/np.sqrt(diag[j])
            work-=v[:,None]*v[None,:];cr-=v[:,None]*h[None,:]
        P=P[piv];u=u[piv]
        w=np.linalg.lstsq(P.T,X.T,rcond=1e-12)[0]
        X=w.T@P
    elif len(u)<a.basis:
        extra=a.basis-len(u);u=np.r_[u,np.eye(n)[np.zeros(extra,dtype=int)]]
        w=np.r_[w,np.zeros((extra,w.shape[1]))]
    val,vec=np.linalg.eigh(X@X.T);order=np.argsort(val)[::-1][:a.columns]
    rot=vec[:,order];w=w@rot;W=of@rot
    if w.shape[1]<a.columns:
        pad=a.columns-w.shape[1];w=np.c_[w,np.zeros((a.basis,pad))];W=np.r_[W,np.zeros(pad)]
    Us.append(u);ws.append(w);Ws.append(W)
    print('group',start,'retained eig norm',float(val[order].sum()),flush=True)
uu=np.array(Us);ww=np.array(ws);WW=np.array(Ws)
np.savez(a.out,U=uu,w=ww,W=WW)
P=uu[:,:,ii]*uu[:,:,jj]*fac;X=np.einsum('rbc,rbp->rcp',ww,P).reshape(-1,len(ii))
g=np.load(HERE/'parameters/chem_reference.npz')['G'];g=g[ix[:,None],ix[None,:]].astype(float)*fac[:,None]*fac[None,:]
print(json.dumps(dict(shape=ww.shape,error_2=float(np.linalg.norm(X.T@X-g)))),flush=True)
