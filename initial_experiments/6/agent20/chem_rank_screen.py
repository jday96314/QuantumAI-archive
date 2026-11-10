"""Screen exact-word rank reduction with a physical generator warm gauge.

This is a physical numerical experiment, not a native submission acceptance.
Reuses read-only chem-audit's native signed target and physical formulas.
"""
import json
import time
from pathlib import Path

import numpy as np
import torch
from audit_core import Metric

HERE=Path(__file__).resolve().parent
BANK=HERE.parent.parent/'5/chem-census/candidates/structured/submission.npz'


def quantities(m,u,w,off):
    with torch.no_grad():
        gap,lam,effective=m.physical(u,w,off)
        return dict(error2=float(m.loss(u,w).sqrt()),gap=float(gap),Lambda=float(lam),
                    lambda_eff=float(effective.sqrt()),weight_norm=float(w.norm()))


if __name__=='__main__':
    torch.set_num_threads(1)
    started=time.time();z=dict(np.load(BANK));m=Metric(z['tree_edges'])
    words=z['angle_words'][76:]
    u=m.decode(m.tensor(words)[None]);w0=m.tensor(z['W'])
    off0=m.tensor(z['H2_I']-z['W'].sum(1))
    reports=[dict(C=350,mode='accepted',**quantities(m,u,w0,off0))]
    print(json.dumps(reports[-1]),flush=True)
    with torch.no_grad():
        P=m.feats(u)[0];Q,R=torch.linalg.qr(P.T)
        H=Q.T@m.G@Q;ev,V=torch.linalg.eigh((H+H.T)/2)
        X0=w0[0].T@P
        left,singular,_=torch.linalg.svd(X0,full_matrices=False)
        physical_offset=off0[0]@X0
        for C in [348,344,342,340]:
            warm=w0[0]@left[:,:C]
            candidate=torch.linalg.solve(R,V[:,-C:]*ev[-C:].clamp(min=0).sqrt())
            a,_,b=torch.linalg.svd(candidate.T@warm,full_matrices=False)
            candidate=candidate@(a@b)
            w=candidate[None]
            X=candidate.T@P
            off=torch.linalg.solve(X@X.T,X@physical_offset)[None]
            before=quantities(m,u,w,off)
            m.project_gap(u,w,off)
            row=dict(C=C,mode='exact_rank_optimum',before=before,
                     after=quantities(m,u,w,off),seconds=time.time()-started)
            reports.append(row);print(json.dumps(row),flush=True)
            np.savez(HERE/'parameters'/f'chem_rank_C{C}.npz',U=u.cpu().numpy(),
                     w=w.cpu().numpy(),W=off.cpu().numpy(),tree_edges=z['tree_edges'],
                     angle_words=words,brot=8)
            (HERE/'results/chem_rank_screen.json').write_text(json.dumps(reports,indent=2)+'\n')
