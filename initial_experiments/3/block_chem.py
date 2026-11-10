"""Block least squares over PSD ensemble matrices with fixed vector rows.

Each block has orthonormal rows Q. Its exact coordinate optimum is the
positive part of Q times the residual-without-that-block times Q.T.
"""
import argparse,time,json
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--cycles',type=int,default=1000);p.add_argument('--columns',type=int)
p.add_argument('--quantize',type=int,default=12)
a=p.parse_args();torch.set_num_threads(3);HERE=Path(__file__).resolve().parent
z=np.load(a.init);u=z['U'].copy();w=z['w'].copy();W=z['W'];R,B,C=w.shape;n=76
base=np.load(HERE/'candidates/chem35c45_rot.npz');tree=base['tree_edges']
if a.quantize:
    from compress_chem import encode
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
    words=encode(u.reshape(-1,n),tree,a.quantize)
    u=reconstruct_u_vectors(words/(1<<(a.quantize-1)),tree).reshape(R,B,n)
else:words=np.empty((0,75),dtype=int)
ii,jj=np.tril_indices(n);fac=np.where(ii==jj,1,np.sqrt(2.));ix=ii*n+jj
P=u[:,:,ii]*u[:,:,jj]*fac
gg=P@P.transpose(0,2,1);ev,vec=np.linalg.eigh(gg)
ev=np.maximum(ev,1e-12)
root=np.einsum('rbi,ri,rci->rbc',vec,np.sqrt(ev),vec)
iroot=np.einsum('rbi,ri,rci->rbc',vec,1/np.sqrt(ev),vec)
Q=iroot@P
h=(root@w)@(root@w).transpose(0,2,1)
g=np.load(HERE/'parameters/chem_reference.npz')['G']
g=g[ix[:,None],ix[None,:]].astype(float)*fac[:,None]*fac[None,:]
d=dict(device='cuda',dtype=torch.float64)
qt=torch.tensor(Q.reshape(R*B,-1),**d);G=torch.tensor(g,**d)
OV=(qt@qt.T).reshape(R,B,R,B).permute(0,2,1,3).contiguous()
T=(qt@G@qt.T).reshape(R,B,R,B)
T=torch.stack([T[r,:,r,:] for r in range(R)])
H=torch.tensor(h,**d);H0=H.clone()
RES=T-torch.einsum('rsab,sbc,rsdc->rad',OV,H,OV)
const=float(G.square().sum());start=time.time();best=float('inf')
for cycle in range(a.cycles):
    for r in range(R):
        target=.5*(RES[r]+RES[r].T)+H[r]
        e,v=torch.linalg.eigh(target)
        e=e.clamp(min=0)
        if a.columns and a.columns<B:e[:-a.columns]=0
        new=(v*e)@v.T
        change=new-H[r];H[r]=new
        A=OV[:,r]
        RES-=A@change@A.transpose(1,2)
    if cycle%10==0 or cycle==a.cycles-1:
        error2=const-float((H*(T+RES)).sum())
        err=np.sqrt(max(0,error2))
        if err<best:
            best=err
            hh=H.cpu().numpy();e,v=np.linalg.eigh(hh);e=np.maximum(e,0)
            Cnew=a.columns or B
            vv=v[:,:,-Cnew:]*np.sqrt(e[:,-Cnew:])[:,None,:]
            ww=iroot@vv
            # Preserve the old identity-offset contribution projected onto
            # the new columns. One-body rows are re-encoded separately.
            oldX=np.einsum('rbc,rbp->rcp',w,P)
            newX=np.einsum('rbc,rbp->rcp',ww,P)
            off=np.einsum('rc,rcp->rp',W,oldX)
            Wnew=np.stack([np.linalg.lstsq(newX[r].T,off[r],rcond=1e-10)[0] for r in range(R)])
            np.savez(a.out,U=u,w=ww,W=Wnew,angle_words=words,tree_edges=tree,err2=err,brot=a.quantize)
        print(json.dumps(dict(cycle=cycle,error_2=err,best=best,seconds=time.time()-start)),flush=True)
