"""Metric audit, orthogonal generator mixing, commutator clustering and JAD.

All X_c are full symmetric orbital matrices: inter-orbital Coulomb entries
are retained. Generator rotations also rotate their scalar identity offsets.
"""
import argparse
import time
import numpy as np
import torch
from scipy.optimize import linear_sum_assignment
from common import HERE,load,feature,unpack,report

def balanced_clusters(distance,R,C,seed):
    """Capacity-constrained medoid/mean-distance assignment, deterministic."""
    n=len(distance);rng=np.random.default_rng(seed)
    centers=[int(np.argmin(distance.sum(0)))]
    while len(centers)<R:
        d=distance[:,centers].min(1);d[centers]=-1
        centers.append(int(d.argmax()))
    slots=np.repeat(np.arange(R),C)[:n]
    cost=distance[:,centers]
    _,col=linear_sum_assignment(cost[:,slots]);labels=slots[col]
    best=float('inf');bestlab=labels.copy()
    for _ in range(30):
        cost=np.stack([distance[:,labels==r].mean(1) for r in range(R)],1)
        # Keep a slight self-attraction to avoid equally good two-cycles.
        cost[np.arange(n),labels]-=1e-9
        _,col=linear_sum_assignment(cost[:,slots]);labels=slots[col]
        val=sum(distance[np.ix_(labels==r,labels==r)].sum() for r in range(R))
        if val<best:best=val;bestlab=labels.copy()
    return [np.flatnonzero(bestlab==r) for r in range(R)]

def commutators(X):
    """||[A,B]||_F^2 / (||A||_F^2 ||B||_F^2), no orbital deletion."""
    norm=X.square().sum((1,2));parts=[]
    for lo in range(0,len(X),8):
        AB=X[lo:lo+8,None]@X[None]
        comm=AB-AB.transpose(-1,-2)
        parts.append(comm.square().sum((-1,-2))/(norm[lo:lo+8,None]*norm[None,:]).clamp(min=1e-30))
    return torch.cat(parts).cpu().numpy()

def jad(X,steps=350):
    """Orthogonal shared basis via differentiable QR; diagonal energy objective.

    Identity shifts have no effect on offdiagonal energy. All generators are
    included, weighted by their physical magnitudes rather than normalized.
    """
    d=dict(device=X.device,dtype=X.dtype)
    e,q=torch.linalg.eigh((X@X).sum(0))
    q=q.clone().requires_grad_()
    opt=torch.optim.Adam([q],lr=.012)
    best=-1.;bestq=q.detach().clone()
    for step in range(steps):
        opt.zero_grad();Q,_=torch.linalg.qr(q)
        diag=torch.einsum('pi,cpq,qi->ci',Q,X,Q)
        objective=-diag.square().sum()
        objective.backward();opt.step()
        energy=float(-objective.detach())
        if energy>best:best=energy;bestq=Q.detach().clone()
    norm=float(X.square().sum())
    return bestq.T,dict(total_factor_norm_squared=norm,offdiagonal_norm_squared=norm-best,relative_offdiagonal=np.sqrt(max(0,norm-best)/norm))

def omp(candidates,target,B,ii,jj,fac):
    """Simultaneous OMP, with exact symmetric projector Gram whitening.

    Candidate pool combines the JAD basis, generator eigenbases and the entire
    incumbent dictionary. Atoms may be nonorthogonal and overcomplete.
    """
    d=dict(device=target.device,dtype=target.dtype)
    cand=candidates/candidates.norm(dim=1,keepdim=True)
    P=cand[:,ii]*cand[:,jj]*fac
    # SOMP tracks orthogonalized candidates without a pool-size-squared Gram.
    rhs=P@target.T;diag=P.square().sum(1);selected=[];Q=[]
    for k in range(B):
        gain=rhs.square().sum(1)/diag.clamp(min=1e-20)
        gain[diag<1e-9]=-1
        if selected:gain[selected]=-1
        j=int(gain.argmax());selected.append(j)
        v=P[j].clone()
        if Q:
            qs=torch.stack(Q);v-=((qs@v)[:,None]*qs).sum(0)
            v-=((qs@v)[:,None]*qs).sum(0)
        v/=v.norm();Q.append(v)
        corr=P@v;rhs-=corr[:,None]*(v@target.T)[None,:]
        diag=(diag-corr.square()).clamp(min=0)
    U=cand[selected]
    Pf=U[:,ii]*U[:,jj]*fac
    gram=Pf@Pf.T
    w=torch.linalg.solve(gram+1e-11*torch.eye(B,**d),Pf@target.T)
    residual=w.T@Pf-target
    return U,w,dict(generator_projection_error=float(residual.norm()),gram_condition=float(torch.linalg.cond(gram)),max_weight=float(w.abs().max()))

def main():
    p=argparse.ArgumentParser();p.add_argument('--R',type=int,default=4);p.add_argument('--B',type=int,default=180);p.add_argument('--C',type=int,default=88)
    p.add_argument('--jad-steps',type=int,default=350);p.add_argument('--tag',default='r4b180c88')
    p.add_argument('--gauge',choices=['commutator_best','physical_PCA','incumbent','random_orthogonal'],default='commutator_best');a=p.parse_args()
    begin=time.time();torch.set_num_threads(2);z,sp=load();d=dict(device='cuda',dtype=torch.float64)
    U=torch.tensor(z['U'][0],**d);w=torch.tensor(z['w'][0],**d);off=torch.tensor(z['offset'][0],**d)
    ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d)
    P=U[:,ii]*U[:,jj]*fac
    F=w.T@P
    G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d)
    X=torch.tensor(unpack(F.cpu().numpy(),sp),**d)
    gram_discrepancy=float(((U@U.T).square()-P@P.T).abs().max())
    packed_norm=float(F.square().sum());matrix_norm=float(X.square().sum())
    assert gram_discrepancy<1e-12 and abs(packed_norm-matrix_norm)<1e-10
    direct=float((F.T@F-G).norm())
    assert abs(direct-.0199953044905532)<1e-8
    e,V=torch.linalg.eigh(F@F.T);V=V.flip(1)
    rng=np.random.default_rng(58431);random=np.linalg.qr(rng.standard_normal((350,350)))[0]
    gauges=[('incumbent',torch.eye(350,**d)),('physical_PCA',V),('random_orthogonal',torch.tensor(random,**d))]
    controls=[];best=None
    for name,O in gauges:
        f=O.T@F;x=torch.einsum('ac,apq->cpq',O,X);o=off@O
        # These checks bind the offset rotation convention to the factors.
        tensor_delta=float((f.T@f-F.T@F).abs().max())
        offset_delta=float((o@f-off@F).abs().max())
        norm_delta=abs(float(o.square().sum()-off.square().sum()))
        assert max(tensor_delta,offset_delta,norm_delta)<1e-9
        dist=commutators(x)
        groups=balanced_clusters(dist,a.R,a.C,519)
        # Use physical-magnitude weighting in selecting the gauge, not just
        # normalized commutators of very small PCA tail columns.
        norms=f.square().sum(1).cpu().numpy()
        weighted=dist*norms[:,None]*norms[None,:]
        within=float(sum(weighted[np.ix_(g,g)].sum() for g in groups))
        entry=dict(gauge=name,tensor_gram_max_discrepancy=tensor_delta,offset_contraction_max_discrepancy=offset_delta,offset_squared_norm_discrepancy=norm_delta,
                   mean_normalized_commutator=float(dist.mean()),within_group_weighted_commutator=within,group_sizes=list(map(len,groups)))
        controls.append(entry)
        print(__import__('json').dumps(entry),flush=True)
        if a.gauge=='commutator_best':
            if best is None or within<best[0]:best=(within,name,O,f,x,o,groups)
        elif name==a.gauge:best=(within,name,O,f,x,o,groups)
    _,name,O,f,x,o,groups=best
    outU=[];outw=[];outo=[];records=[]
    for r,inds in enumerate(groups):
        xr=x[inds];target=f[inds]
        basis,jadrec=jad(xr,a.jad_steps)
        # JAD alone is screened; overcomplete dictionary fitting proceeds if
        # it cannot retain the full correlated generators.
        _,eigvec=torch.linalg.eigh(xr)
        candidates=torch.cat([basis,U,eigvec.transpose(1,2).reshape(-1,76)])
        uu,ww,omprec=omp(candidates,target,a.B,ii,jj,fac)
        pad=a.C-len(inds)
        ww=torch.nn.functional.pad(ww,(0,pad))
        oo=torch.nn.functional.pad(o[inds],(0,pad))
        outU.append(uu);outw.append(ww);outo.append(oo)
        records.append(dict(r=r,columns=inds.tolist(),jad=jadrec,overcomplete=omprec))
        print(__import__('json').dumps(dict(r=r,jad=jadrec,overcomplete=omprec)),flush=True)
    un=torch.stack(outU);wn=torch.stack(outw);on=torch.stack(outo)
    fn=(wn.transpose(1,2)@(un[:,:,ii]*un[:,:,jj]*fac)).reshape(a.R*a.C,-1)
    err=float((fn.T@fn-G).norm())
    np.savez(HERE/f'parameters/{a.tag}_seed.npz',U=un.cpu().numpy(),w=wn.cpu().numpy(),W=on.cpu().numpy(),tree_edges=z['tree_edges'],err2=err)
    np.savez(HERE/f'parameters/{a.tag}_generator_gauge.npz',orthogonal=O.cpu().numpy(),offsets=o.cpu().numpy(),generator_features=f.cpu().numpy())
    report(HERE/f'results/{a.tag}_generators.json',dict(status='CONTINUOUS_INITIALIZER_ONLY',R=a.R,B=a.B,C=a.C,
        incumbent_error_2=direct,projector_gram_max_discrepancy=gram_discrepancy,packed_factor_norm_squared=packed_norm,matrix_factor_norm_squared=matrix_norm,
        gauges=controls,chosen_gauge=name,ensembles=records,seed_error_2=err,seconds=time.time()-begin))

if __name__=='__main__':main()
