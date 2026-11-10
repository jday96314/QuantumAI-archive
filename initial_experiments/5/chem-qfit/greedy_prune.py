"""Sequential tensor-aware pruning, recomputing leverage after each deletion."""
import argparse,json,time
from core import *
p=argparse.ArgumentParser();p.add_argument('--init',default=str(INC/'submission.npz'));p.add_argument('--B',type=int,default=671)
p.add_argument('--prefix',default='greedy');p.add_argument('--drop',type=int,default=1);a=p.parse_args()
z,u,w,tree,words=load(a.init);m=Metric(tree);ut=m.tensor(u);wt=m.tensor(w);start=time.time()
while len(ut[0])>a.B:
    P=m.feats(ut)[0];S=P@P.T;inv=torch.linalg.inv(S);K=wt[0]@wt[0].T
    qd=inv.diag();score=2*(K@S@K).diag()/qd-(K.diag()/qd).square()
    B=max(a.B,len(ut[0])-a.drop);ix=score.argsort(descending=True)[:B].sort().values
    ut=ut[:,ix];wt=wt[:,ix];words=words[ix.cpu().numpy()];wt,err=m.project(ut,wt,1)
    if B%8==0 or B in [711,700,671,639,600,575,543,511] or B==a.B:
        save(HERE/f'candidates/{a.prefix}_B{B}.npz',ut,wt,tree,words,err2=err)
        print(json.dumps(dict(B=B,err=err,seconds=time.time()-start)),flush=True)
