"""Make integer controls and small single-block dictionaries from selected words."""
import json
from core import *
z,u,w,tree,words=load(INC/'submission.npz');m=Metric(tree)
save(HERE/'candidates/incumbent_control.npz',u,w,tree,words)
ut=m.tensor(u);wt=m.tensor(w)
wt,err=m.project(ut,wt,1)
save(HERE/'candidates/incumbent_projected.npz',ut,wt,tree,words,err2=err)
P=m.feats(ut)[0];S=P@P.T;K=wt[0]@wt[0].T
# Exact leave-one-out projection score for the full rank dictionary. This is
# only a pruning proposal; every resulting rank-C residual is evaluated.
inv=torch.linalg.inv(S);qd=inv.diag();Y=K@S@K
importance=2*Y.diag()/qd-(K.diag()/qd).square()
order=torch.argsort(importance,descending=True)
for B in [700,671,639,600,575,543,511]:
    ix=order[:B].sort().values;uu=ut[:,ix];ww=wt[:,ix];ww,err=m.project(uu,ww,1)
    save(HERE/f'candidates/pruned_B{B}.npz',uu,ww,tree,words[ix.cpu().numpy()],err2=err)
    print(json.dumps(dict(B=B,C=ww.shape[-1],integer_error=err)),flush=True)
