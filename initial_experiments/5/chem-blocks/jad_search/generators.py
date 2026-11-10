"""Derive physical matrices, mix their columns, cluster commutators, JAD, OMP.

No orbital blocks are discarded. Each ensemble retains full 76x76 matrices.
Orthogonal column mixing also rotates native identity offsets.
"""
import argparse,time
import numpy as np
import torch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from core import *

def commutator_dist(mats):
    # Eight fixed random probes estimate ||[X_i,X_j]||_F^2. Identity parts
    # commute and are removed before relative normalization.
    d=dict(device='cuda',dtype=torch.float64)
    M=torch.tensor(mats,**d)
    M=M-torch.diagonal(M,dim1=1,dim2=2).sum(1)[:,None,None]/76*torch.eye(76,**d)
    norms=M.square().sum((1,2)).sqrt().clamp(min=1e-15)
    M=M/norms[:,None,None]
    Z=torch.tensor(np.random.default_rng(531).choice([-1.,1.],size=(76,8)),**d)
    Y=M@Z;rows=[]
    for start in range(0,len(M),20):
        a=torch.einsum('ipq,jqz->ijpz',M[start:start+20],Y)
        b=torch.einsum('jpq,iqz->ijpz',M,Y[start:start+20])
        rows.append(((a-b).square().sum((2,3))/8).cpu().numpy())
    dist=np.concatenate(rows);dist=(dist+dist.T)/2;np.fill_diagonal(dist,0)
    return dist

def balanced_cluster(dist,R,C):
    groups=[[i] for i in range(len(dist))]
    sums=dist.copy();sizes=np.ones(len(dist),dtype=int)
    while len(groups)>R:
        score=sums/(sizes[:,None]*sizes[None,:])
        score[sizes[:,None]+sizes[None,:]>C]=np.inf
        score[np.triu_indices(len(groups))]=np.inf
        i,j=np.unravel_index(np.argmin(score),score.shape)
        if not np.isfinite(score[i,j]):
            # Balanced capacity can trap agglomeration. Pack its ordered leaves
            # into C-sized groups, then improve using capacity-preserving swaps.
            order=sum(sorted(groups,key=lambda g:g[0]),[])
            groups=[order[i*C:(i+1)*C] for i in range(R)]
            break
        groups[j]+=groups[i];groups.pop(i)
        sums[j,:]+=sums[i,:];sums[:,j]+=sums[:,i];sizes[j]+=sizes[i]
        sums=np.delete(np.delete(sums,i,axis=0),i,axis=1);sizes=np.delete(sizes,i)
    def objective(gs):return sum(dist[np.ix_(g,g)].sum() for g in gs)
    rng=np.random.default_rng(792);best=objective(groups)
    for _ in range(12000):
        r,s=rng.choice(R,2,replace=False);i=int(rng.integers(len(groups[r])));j=int(rng.integers(len(groups[s])))
        u,v=groups[r][i],groups[s][j]
        delta=2*(dist[v,groups[r]].sum()-dist[u,groups[r]].sum()+dist[u,groups[s]].sum()-dist[v,groups[s]].sum()-2*dist[u,v])
        if delta<-1e-14:
            groups[r][i],groups[s][j]=v,u;best+=delta
    return groups

def jad(matrices,sweeps):
    A=matrices.copy();Q=np.eye(76);records=[]
    initial=float((A*A).sum()-np.diagonal(A,axis1=1,axis2=2).__pow__(2).sum())
    for sweep in range(sweeps):
      for p in range(76):
       for q in range(p):
        a=A[:,p,p]-A[:,q,q];b=2*A[:,p,q]
        aa=a@a;bb=b@b;ab=a@b
        theta=.25*np.arctan2(2*ab,aa-bb)
        if abs(theta)<1e-10:continue
        c,s=np.cos(theta),np.sin(theta)
        # Q <- Q R, R=[[c,-s],[s,c]], A <- R.T A R.
        ap,aq=A[:,:,p].copy(),A[:,:,q].copy();A[:,:,p]=c*ap+s*aq;A[:,:,q]=-s*ap+c*aq
        ap,aq=A[:,p,:].copy(),A[:,q,:].copy();A[:,p,:]=c*ap+s*aq;A[:,q,:]=-s*ap+c*aq
        qp,qq=Q[:,p].copy(),Q[:,q].copy();Q[:,p]=c*qp+s*qq;Q[:,q]=-s*qp+c*qq
      residual=float((A*A).sum()-np.diagonal(A,axis1=1,axis2=2).__pow__(2).sum())
      records.append(residual)
    assert np.max(np.abs(Q.T@Q-np.eye(76)))<1e-12
    return Q.T,dict(initial_offdiag_squared=initial,final_offdiag_squared=residual,sweeps=records)

def omp(candidates,target,B):
    d=dict(device='cuda',dtype=torch.float64)
    U=torch.tensor(candidates,**d)
    sp=dict(np.load(SPECTRUM));P=U[:,sp['ii']]*U[:,sp['jj']]*torch.tensor(sp['fac'],**d)
    T=torch.tensor(target,**d)
    K=P@P.T;K0=K.clone();rhs=P@T.T;selected=[]
    for _ in range(B):
        diag=K.diag().clamp(min=1e-16)
        gain=rhs.square().sum(1)/diag;gain[diag<1e-9]=-1
        if selected:gain[selected]=-1
        j=int(gain.argmax());selected.append(j)
        v=K[:,j]/diag[j].sqrt();h=rhs[j]/diag[j].sqrt()
        K-=v[:,None]*v[None,:];rhs-=v[:,None]*h[None,:]
    ix=torch.tensor(selected,device='cuda')
    w=torch.linalg.solve(K0[ix][:,ix]+1e-12*torch.eye(B,**d),P[ix]@T.T)
    return candidates[selected],w.cpu().numpy(),float((w.T@P[ix]-T).square().sum().sqrt())

def main():
    p=argparse.ArgumentParser();p.add_argument('--R',type=int,default=4);p.add_argument('--B',type=int,default=180);p.add_argument('--C',type=int,default=88)
    p.add_argument('--sweeps',type=int,default=6);p.add_argument('--mix',choices=['bestcomm','original','spectral','random'],default='bestcomm');a=p.parse_args();torch.set_num_threads(2);beg=time.time()
    z,sp=load();P=feature(z['U'][0],sp);X=z['w'][0].T@P;off=z['offset'].ravel();G=gram(sp)
    M=unpack(X,sp);K=P@P.T
    checks=dict(projector_gram_error=float(np.max(np.abs(K-(z['U'][0]@z['U'][0].T)**2))),
        matrix_feature_norm_error=float(abs((M*M).sum()-(X*X).sum())),
        incumbent_error_2=float(np.linalg.norm(X.T@X-G)),
        generator_reconstruction_error=float(np.max(np.abs(M-np.einsum('bc,bp,bq->cpq',z['w'][0],z['U'][0],z['U'][0])))))
    assert checks['projector_gram_error']<1e-12 and checks['generator_reconstruction_error']<1e-12
    # Metric-aware covariance eigensystem: K=L L^T, H=L^T w w^T L.
    L=np.linalg.cholesky(K+1e-12*np.eye(len(K)));Qp=np.linalg.solve(L,P)
    H=L.T@z['w'][0]@z['w'][0].T@L
    checks['whitened_covariance_error']=float(np.linalg.norm(Qp.T@H@Qp-X.T@X))
    modes={}
    e,V=np.linalg.eigh(X@X.T)
    modes['spectral']=V[:,::-1].T
    modes['original']=np.eye(len(X))
    rng=np.random.default_rng(410);q,_=np.linalg.qr(rng.standard_normal((len(X),len(X))))
    modes['random']=q.T
    best=None;controls=[]
    for name,O in modes.items():
        if a.mix!='bestcomm' and name!=a.mix:continue
        Y=O@X;shift=O@off
        check=dict(mode=name,gram_invariance=float(np.linalg.norm(Y.T@Y-X.T@X)),offset_invariance=float(np.linalg.norm(shift@Y-off@X)),offset_norm_error=float(abs(shift@shift-off@off)))
        dist=commutator_dist(unpack(Y,sp));groups=balanced_cluster(dist,a.R,a.C)
        check['commutator_within_mean']=float(np.mean([dist[np.ix_(g,g)].mean() for g in groups]));check['commutator_all_mean']=float(dist.mean())
        controls.append(check);emit(check)
        if best is None or check['commutator_within_mean']<best[0]:best=(check['commutator_within_mean'],name,Y,shift,groups)
    _,name,Y,shift,groups=best
    outs=[];weights=[];offsets=[];jadrows=[];targets=[]
    for r,group in enumerate(groups):
        target=Y[group];Uj,jr=jad(unpack(target,sp),a.sweeps)
        # Overcomplete nonorthogonal candidates include full-orbital incumbent
        # atoms and eigenvectors of several high-energy generators.
        _,vec=np.linalg.eigh(unpack(target[np.argsort(np.linalg.norm(target,axis=1))[-8:]],sp))
        candidates=np.concatenate([Uj,z['U'][0],vec.transpose(0,2,1).reshape(-1,76)])
        U,w,local=omp(candidates,target,a.B)
        pad=a.C-len(group);w=np.pad(w,((0,0),(0,pad)))
        offs=np.pad(shift[group],(0,pad));t=np.pad(target,((0,pad),(0,0)))
        outs.append(U);weights.append(w);offsets.append(offs);targets.append(t)
        jr.update(group=r,C_actual=len(group),factor_projection_error=local);jadrows.append(jr);emit(jr)
    U=np.array(outs);w=np.array(weights);offset=np.array(offsets);T=np.array(targets)
    model=np.einsum('rbc,rbd->rcd',w,feature(U,sp)).reshape(a.R*a.C,-1)
    err=float(np.linalg.norm(model.T@model-G))
    stem=f'r{a.R}_b{a.B}_c{a.C}'+('' if a.mix=='bestcomm' else '_'+a.mix)
    np.savez(HERE/'candidates'/f'{stem}_seed.npz',U=U,w=w,offset=offset,target=T,tree_edges=z['tree_edges'],err2=err)
    np.savez(HERE/'candidates'/'incumbent_generators.npz',X=X,matrices=M,offset=off,feature_gram=K)
    data=dict(status='CONTINUOUS_INITIALIZER_ONLY',checks=checks,mixing_controls=controls,selected_mixing=name,groups=groups,jad=jadrows,error_2=err,shape=[a.R,a.B,a.C],seconds=time.time()-beg)
    report(f'{stem}_generators.json',data);emit(data)

if __name__=='__main__':main()
