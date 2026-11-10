"""Capacity screen for exact low-bit sparse orbital-projector dictionaries.

A block contains all76 coordinate projectors plus K pair projectors
(e_i+e_j)/sqrt2. Its feature span is exactly the76 diagonal coordinates plus
those K off-diagonal coordinates. Three-bit Givens words suffice exactly on
any binary aggregation tree. For disjoint pair groups, every cross-group
pair/pair entry of the signed target is forced to zero, yielding a rigorous
fixed-support residual lower bound before PSD/rank/normalization restrictions.
This screen changes neither reference data nor host acceptance.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
from audit_core import Metric
from chem_rank_screen import BANK
HERE=Path(__file__).resolve().parent


def greedy(A,capacity,seed):
    rng=np.random.default_rng(seed);n=len(A);unused=np.ones(n,dtype=bool);groups=[]
    degree=A.sum(1)
    while unused.any():
        if seed==0:j=int(np.where(unused,degree,-1).argmax())
        else:
            candidates=np.flatnonzero(unused)
            weights=degree[candidates]+1e-30
            j=int(rng.choice(candidates,p=weights/weights.sum()))
        unused[j]=False;group=[j];score=A[:,j].copy()
        for _ in range(min(capacity-1,int(unused.sum()))):
            j=int(np.where(unused,score,-1).argmax());unused[j]=False;group.append(j);score+=A[:,j]
        groups.append(np.asarray(group,dtype=int))
    return groups


def main():
    p=argparse.ArgumentParser();p.add_argument('--capacity',type=int,nargs='+',default=[52,104,152,228,304]);p.add_argument('--seeds',type=int,default=8);a=p.parse_args()
    started=time.time();z=dict(np.load(BANK));m=Metric(z['tree_edges'],device='cpu')
    G=m.G.numpy();ii=m.ii.numpy();jj=m.jj.numpy();off=np.flatnonzero(ii!=jj);diag=np.flatnonzero(ii==jj)
    Go=G[np.ix_(off,off)];A=Go.square() if hasattr(Go,'square') else Go**2
    total=float(A.sum());np.fill_diagonal(A,0);rows=[]
    for K in a.capacity:
        for seed in range(a.seeds):
            groups=greedy(A,K,seed)
            inside=sum(float((Go[np.ix_(g,g)]**2).sum()) for g in groups)
            bound=np.sqrt(max(0.,total-inside));R=len(groups)
            row=dict(K=K,B=76+K,R=R,RC_lower_rank=290,C_min_from_rank=int(np.ceil(290/R)),
                     seed=seed,forced_error_lower_bound=float(bound),seconds=time.time()-started)
            rows.append(row);print(json.dumps(row),flush=True)
            memberships=np.full(len(off),-1,dtype=int)
            for r,g in enumerate(groups):memberships[g]=r
            np.savez(HERE/'parameters'/f'sparse_pairs_K{K}_s{seed}.npz',offdiag_pair_indices=off,
                     group=memberships,ii=ii,jj=jj,K=K)
        (HERE/'results/sparse_pair_capacity_screen.json').write_text(json.dumps(dict(status='FINITE_FIXED_SUPPORT_SCREEN',arguments=vars(a),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            offdiag_tensor_norm=float(np.linalg.norm(Go)),all_target_norm=float(np.linalg.norm(G)),rows=rows,
            scope='rigorous forced cross-block residual for each saved disjoint support family; necessary bound only; not a global sparse/overlapping-dictionary bound'),indent=2)+'\n')

if __name__=='__main__':main()
