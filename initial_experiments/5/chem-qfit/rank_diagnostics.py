"""Separate fixed integer dictionary capacity, rank C, and block solver error."""
import json
from core import *
for tag in ['qd_b479','qd_b447']:
    z,u,w,tree,words=load(ROOT/f'initial_experiments/4/claude/chem/{tag}.npz');m=Metric(tree);ut=m.tensor(u);wt=m.tensor(w)
    P=m.feats(ut).flatten(0,1);Q,_=torch.linalg.qr(P.T);H=Q.T@m.G@Q
    full_span=float((m.G-Q@H@Q.T).norm())
    for C in [160,175,180,190,200,225]:
        ww=torch.zeros((2,w.shape[1],C),**m.d);ww[:,:,:min(C,175)]=wt[:,:,:min(C,175)]
        ww,err=m.project(ut,ww,12)
        save(HERE/f'candidates/{tag}_C{C}.npz',ut,ww,tree,words,err2=err)
        print(json.dumps(dict(tag=tag,C=C,full_joint_span_lower_bound=full_span,block_error=err)),flush=True)
