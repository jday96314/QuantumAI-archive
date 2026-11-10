"""Non-JAD initializations: accepted covariance graph and residual eigenatoms.

Graph construction shares high-coupling accepted atoms across ensembles and
balances exclusive capacity. Eigenatom construction uses no common eigenbasis:
each target generator contributes its individual full-orbital eigenvectors;
simultaneous pivoting selects a nonorthogonal projector span for each group.
"""
import argparse
from scipy.optimize import linear_sum_assignment
from scipy.linalg import eigh
from core import *

def balance(embed,R):
    n=len(embed);centers=[int(np.argmax((embed*embed).sum(1)))]
    for _ in range(1,R):
        dist=((embed[:,None]-embed[centers][None])**2).sum(2).min(1);dist[centers]=-1
        centers.append(int(dist.argmax()))
    slots=np.arange(n)%R;c=embed[centers]
    labels=None
    for _ in range(30):
        d=((embed[:,None]-c[None])**2).sum(2)
        _,assign=linear_sum_assignment(d[:,slots]);new=slots[assign]
        if labels is not None and np.array_equal(labels,new):break
        labels=new;c=np.stack([embed[labels==r].mean(0) for r in range(R)])
    return labels

def graph(R,B,shared):
    z=dict(np.load(INC));u=reconstruct_u_vectors(z['angle_words'][76:]/128,z['tree_edges']);w=z['W'][0]
    K=w@w.T;G=(u@u.T)**2
    # Atom attribution includes covariance coupling, rather than raw row size.
    score=np.abs(K*G).sum(1)
    count=min(len(u),R*(B-shared)+shared)
    chosen=np.argsort(score)[-count:]
    coupling=K[np.ix_(chosen,chosen)]**2
    degree=coupling.sum(1);norm=coupling/np.sqrt(degree[:,None]*degree[None,:]+1e-30)
    _,V=eigh(norm,subset_by_index=[len(chosen)-R,len(chosen)-1]);labels=balance(V,R)
    cross=np.stack([coupling[:,labels==r].sum(1) for r in range(R)],1)
    bridge=(cross.sum(1)-cross.max(1))*np.sqrt(score[chosen]+1e-30)
    common=chosen[np.argsort(bridge)[-shared:]] if shared else np.array([],dtype=int)
    rem=np.setdiff1d(chosen,common);A=K[np.ix_(rem,rem)]**2;degree=A.sum(1)
    _,V=eigh(A/np.sqrt(degree[:,None]*degree[None,:]+1e-30),subset_by_index=[len(rem)-R,len(rem)-1]);lab=balance(V,R)
    groups=[]
    for r in range(R):
        inds=list(common)+list(rem[lab==r])
        # Fill a smaller exclusive group with the strongest remaining atoms.
        for j in np.argsort(score)[::-1]:
            if len(inds)>=B:break
            if j not in inds:inds.append(int(j))
        groups.append([int(i) for i in inds[:B]])
    return u[np.array(groups)],dict(strategy='accepted_covariance_graph',shared=shared,atom_indices=groups,union_atoms=len(set(sum(groups,[]))))

@torch.no_grad()
def pivot(pool,F,B,common=None):
    P=tensor(features(pool));F=tensor(F);rhs=P@F.T;diag=P.square().sum(1);Q=[];selected=[]
    for k in range(B):
        gain=rhs.square().sum(1)/diag.clamp(min=1e-20);gain[diag<1e-10]=-1
        if selected:gain[selected]=-1
        j=common[k] if common is not None and k<len(common) else int(gain.argmax())
        if float(diag[j])<1e-10:j=int(gain.argmax())
        selected.append(j);v=P[j].clone()
        if Q:
            qq=torch.stack(Q);v-=qq.T@(qq@v);v-=qq.T@(qq@v)
        v/=v.norm();Q.append(v);corr=P@v;rhs-=corr[:,None]*(F@v)[None,:];diag=(diag-corr*corr).clamp(min=0)
    return pool[selected],selected

def eigenatoms(R,B,shared):
    ev,V=eigh(target(),subset_by_index=[2926-360,2925]);order=np.argsort(ev)[::-1];ev,V=ev[order],V[:,order]
    F=V.T*np.sqrt(ev.clip(0))[:,None]
    matrices=unpack(F);_,vectors=np.linalg.eigh(matrices[:180])
    z=dict(np.load(INC));accepted=reconstruct_u_vectors(z['angle_words'][76:]/128,z['tree_edges'])
    pool=np.r_[accepted,vectors.transpose(0,2,1).reshape(-1,76)]
    # Spatial participation provides a materially different grouping from
    # commutator/JAD clustering; every matrix retains all orbital entries.
    orbital_energy=(matrices*matrices).sum(2);labels=balance(orbital_energy/(orbital_energy.sum(1,keepdims=True)+1e-30),R)
    common=[]
    if shared:
        _,common=pivot(pool,F,shared)
    out=[];records=[]
    for r in range(R):
        u,inds=pivot(pool,F[labels==r],B,common);out.append(u);records.append(inds)
    return np.stack(out),dict(strategy='individual_target_eigenatoms_spatial_grouping',shared=shared,pool_size=len(pool),pool_indices=records)

def main():
    p=argparse.ArgumentParser();p.add_argument('--R',type=int,required=True);p.add_argument('--B',type=int,required=True)
    p.add_argument('--shared',type=int,default=32);p.add_argument('--strategy',choices=['graph','eigenatoms'],required=True);p.add_argument('--tag',required=True);a=p.parse_args()
    u,info=(graph if a.strategy=='graph' else eigenatoms)(a.R,a.B,a.shared)
    np.savez(HERE/f'parameters/{a.tag}.npz',U=u,tree_edges=np.load(INC)['tree_edges'])
    report(HERE/f'results/{a.tag}_construction.json',dict(R=a.R,B=a.B,**info))
if __name__=='__main__':main()
