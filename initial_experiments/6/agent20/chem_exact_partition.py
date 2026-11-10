"""Exact-word multiblock bootstrap by Ky Fan column rotations and joint OMP.

Uses the immutable accepted dictionary; no continuous angle quantization.
For two fixed spans, the best balanced generator-column partition is the
top/bottom eigenspace of X(P0-P1)X^T.  Greedy atom updates use the exact
projector Gram metric.  This is a numerical screen, not native acceptance.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric
from chem_rank_screen import BANK,quantities

HERE=Path(__file__).resolve().parent

@torch.no_grad()
def joint_omp(S,correlation,B):
    gram=S.clone();cross=correlation.clone();chosen=[]
    used=torch.zeros(len(S),device=S.device,dtype=torch.bool)
    for _ in range(B):
        diag=gram.diag().clamp(min=1e-13)
        scores=cross.square().sum(1)/diag
        scores[used]=-1
        j=int(scores.argmax());chosen.append(j);used[j]=True
        col=gram[:,j].clone()/diag[j].sqrt()
        target=cross[j].clone()/diag[j].sqrt()
        gram-=col[:,None]*col[None,:]
        cross-=col[:,None]*target[None,:]
    return torch.tensor(chosen,device=S.device)

def main():
    p=argparse.ArgumentParser();p.add_argument('--B',type=int,default=351)
    p.add_argument('--cycles',type=int,default=5);p.add_argument('--seeds',type=int,default=3)
    p.add_argument('--refit-cycles',type=int,default=10);a=p.parse_args()
    torch.set_num_threads(1);z=dict(np.load(BANK));m=Metric(z['tree_edges'])
    words=z['angle_words'][76:];u=m.decode(m.tensor(words)[None]);P=m.feats(u)[0]
    w=m.tensor(z['W'])[0];off=m.tensor(z['H2_I']-z['W'].sum(1))[0]
    S=P@P.T;X=w.T@P;start=time.time();rows=[]
    for seed in range(a.seeds):
        torch.manual_seed(seed)
        if seed==0:
            # Native generator gauge; physically benign starting point.
            V=torch.eye(350,**m.d)
        else:V=torch.linalg.qr(torch.randn((350,350),**m.d))[0]
        for cycle in range(a.cycles):
            corr=S@w@V
            ids=[joint_omp(S,corr[:,r*175:(r+1)*175],a.B) for r in range(2)]
            qs=[torch.linalg.qr(P[ix].T)[0] for ix in ids]
            scores=[X@q for q in qs]
            H=scores[0]@scores[0].T-scores[1]@scores[1].T
            ev,E=torch.linalg.eigh((H+H.T)/2)
            V=torch.cat([E[:,-175:],E[:,:175]],1)
            uu=torch.stack([u[0,ix] for ix in ids])
            aw=np.stack([words[ix.cpu().numpy()] for ix in ids])
            ww=torch.stack([torch.linalg.lstsq(P[ix].T,(V[:,r*175:(r+1)*175].T@X).T).solution
                            for r,ix in enumerate(ids)])
            oo=(off@V).reshape(2,175)
            # Optimal column rotation has a checkable fixed-span projection loss.
            projection_loss=sum(float((V[:,r*175:(r+1)*175].T@X-
                ww[r].T@P[ids[r]]).square().sum()) for r in range(2))
            before=quantities(m,uu,ww,oo)
            ww,error=m.project(uu,ww,a.refit_cycles)
            m.project_gap(uu,ww,oo)
            after=quantities(m,uu,ww,oo)
            row=dict(seed=seed,cycle=cycle,B=a.B,projection_loss=projection_loss,
                     before=before,after=after,seconds=time.time()-start)
            rows.append(row);print(json.dumps(row),flush=True)
            tag=f'exact_partition_B{a.B}_s{seed}_c{cycle}'
            np.savez(HERE/'parameters'/f'{tag}.npz',U=uu.cpu().numpy(),w=ww.cpu().numpy(),
                     W=oo.cpu().numpy(),tree_edges=z['tree_edges'],angle_words=aw.reshape(-1,75),brot=8)
            (HERE/'results'/f'exact_partition_B{a.B}.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':
    with torch.no_grad():main()
