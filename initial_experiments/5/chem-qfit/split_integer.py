"""Split accepted integer dictionary columns, prune each block without rounding.

Keeps the accepted angle basin throughout; differs from splitting a float
dictionary and quantizing the resulting optimized vectors from scratch.
"""
import argparse,json,time
from scipy.cluster.hierarchy import linkage,leaves_list
from core import *
p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--order',choices=['cluster','alternating','spectral'],default='cluster')
a=p.parse_args();z,u,w,tree,words=load(INC/'submission.npz');m=Metric(tree);C=w.shape[-1]
if a.order=='cluster':
    support=np.abs(w[0]);support/=np.linalg.norm(support,axis=0,keepdims=True)
    ids=leaves_list(linkage(support.T,method='ward'))
elif a.order=='alternating':ids=np.r_[np.arange(0,C,2),np.arange(1,C,2)]
else:
    cov=np.abs(w[0]).T@np.abs(w[0]);d=cov.sum(1);cov/=np.sqrt(d[:,None]*d[None,:]);e,v=np.linalg.eigh(cov);ids=np.argsort(v[:,-2])
groups=[ids[:175],ids[175:]];ut=m.tensor(np.repeat(u,2,axis=0));wt=m.tensor(np.stack([w[0][:,g] for g in groups]));aw=np.stack([words,words]);start=time.time()
for B in [700,671,639,607,575,543,511,495,479,463,447]:
    P=m.feats(ut);keep=[]
    for r in range(2):
        S=P[r]@P[r].T;iv=torch.linalg.inv(S);K=wt[r]@wt[r].T
        score=2*(K@S@K).diag()/iv.diag()-(K.diag()/iv.diag()).square()
        keep.append(score.argsort(descending=True)[:B].sort().values)
    ut=torch.stack([ut[r,keep[r]] for r in range(2)]);wt=torch.stack([wt[r,keep[r]] for r in range(2)])
    aw=np.stack([aw[r,keep[r].cpu().numpy()] for r in range(2)]);wt,err=m.project(ut,wt,6)
    tag=str(Path(a.out).with_suffix(''))+f'_B{B}.npz';save(tag,ut,wt,tree,aw.reshape(-1,75),err2=err)
    print(json.dumps(dict(B=B,error=err,s=time.time()-start)),flush=True)
