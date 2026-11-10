"""One resource-favorable pivot: expand a warm multiblock dictionary.

Uses full-orbital symmetric projectors, with no discarded inter-orbital terms.
Additional directions greedily span a Procrustes-aligned reference square root.
The old dictionary and weights remain an exact feasible initialization.
"""
import argparse,time
import numpy as np
import torch
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--B',type=int,default=224);p.add_argument('--output',required=True);a=p.parse_args()
    torch.set_num_threads(2);beg=time.time();z=dict(np.load(a.input));inc,sp=load();U=z['U'];w=z['w'];R,B,N=U.shape;C=w.shape[2];K=R*C
    d=dict(device='cuda',dtype=torch.float64);fac=torch.tensor(sp['fac'],**d)
    X=torch.tensor(np.einsum('rbc,rbd->rcd',w,feature(U,sp)).reshape(K,-1),**d)
    F=torch.tensor((sp['eigenvectors'][:,:K]*np.sqrt(np.maximum(sp['eigenvalues'][:K],0))).T,**d)
    l,_,rh=torch.linalg.svd(X@F.T,full_matrices=False);Y=(l@rh@F).reshape(R,C,-1)
    outs=[];weights=[];records=[]
    for r in range(R):
        target=Y[r];_,v=np.linalg.eigh(unpack(target.cpu().numpy()[np.argsort(np.linalg.norm(target.cpu().numpy(),axis=1))[-8:]],sp))
        candidates=np.concatenate([U[r],inc['U'][0],v.transpose(0,2,1).reshape(-1,N)])
        P=torch.tensor(feature(candidates,sp),**d);S=P@P.T;S0=S.clone();rhs=P@target.T
        selected=[]
        for j in range(a.B):
            diag=S.diag().clamp(min=1e-16)
            gain=rhs.square().sum(1)/diag;gain[diag<1e-9]=-1
            if selected:gain[selected]=-1
            idx=j if j<B else int(gain.argmax());selected.append(idx)
            vv=S[:,idx]/diag[idx].sqrt();h=rhs[idx]/diag[idx].sqrt();S-=vv[:,None]*vv[None,:];rhs-=vv[:,None]*h[None,:]
        ix=torch.tensor(selected,device='cuda');ww=torch.linalg.solve(S0[ix][:,ix]+1e-12*torch.eye(a.B,**d),P[ix]@target.T)
        outs.append(candidates[selected]);weights.append(ww.cpu().numpy());records.append(dict(r=r,added_indices=selected[B:]))
    un=np.array(outs);wn=np.array(weights);G=gram(sp)
    Xnew=np.einsum('rbc,rbd->rcd',wn,feature(un,sp)).reshape(K,-1);error=float(np.linalg.norm(Xnew.T@Xnew-G))
    old=float(np.linalg.norm(X.cpu().numpy().T@X.cpu().numpy()-G));selected='projected'
    if old<error:wn=np.pad(w,((0,0),(0,a.B-B),(0,0)));error=old;selected='preserved_weights'
    np.savez(a.output,U=un,w=wn,offset=z['offset'],target=Y.cpu().numpy(),tree_edges=z['tree_edges'],err2=error)
    data=dict(status='CONTINUOUS_PIVOT_INITIALIZER',input=a.input,shape=[R,a.B,C],old_error_2=old,projected_error_2=float(np.linalg.norm(Xnew.T@Xnew-G)),selected_initialization=selected,error_2=error,records=records,seconds=time.time()-beg)
    report(Path(a.output).stem+'_grow.json',data);emit(data)

if __name__=='__main__':main()
