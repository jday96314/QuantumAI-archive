"""Matched controls: cyclic PSD projection versus joint whitened factor search.

Words never change. This diagnoses numerical block/rank basins before spending
dictionary proposals against them. Trajectories retain superior intermediates.
"""
import argparse,json,time
from compact import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=400);p.add_argument('--seed',type=int,default=0)
p.add_argument('--mode',choices=['lbfgs','cycles','restart'],default='lbfgs');a=p.parse_args()
torch.manual_seed(a.seed);z,u,w,tree,words=load(a.init);m=Metric(tree);t=m.tensor(words.reshape(*u.shape[:2],75));u=m.decode(t);wt=m.tensor(w)
c=Compact(m,u);f=c.factors(wt);start=time.time();best=float(m.loss(u,wt).sqrt());out=Path(a.out)
save(out,u,wt,tree,t,err2=best)
if a.mode=='restart':
    # Rebalance the two factors inside their actual block spans. Cross-block
    # orthogonal column rotations are followed by block-span projection.
    xf=torch.bmm(c.A,f).transpose(1,2).flatten(0,1).T
    rand=torch.randn((xf.shape[1],xf.shape[1]),**m.d);rot=torch.linalg.qr(torch.eye(len(rand),**m.d)+.1*rand/len(rand)**.5).Q
    xf=xf@rot;f=torch.stack([c.A[r].T@xf[:,r*f.shape[-1]:(r+1)*f.shape[-1]] for r in range(c.R)])
f=f.detach().requires_grad_(a.mode!='cycles')
opt=torch.optim.LBFGS([f],lr=1.,max_iter=8,history_size=25,tolerance_grad=1e-12,tolerance_change=1e-16,line_search_fn='strong_wolfe') if a.mode!='cycles' else None
trajectory=[]
for step in range(a.steps):
    if opt:
        def closure():
            opt.zero_grad();loss=c.loss(f);loss.backward();return loss
        opt.step(closure)
    else:f,er=c.project(f,1)
    if step%10==0 or step==a.steps-1:
        er=float(c.loss(f).sqrt());ww=c.weights(f).detach();exact=float(m.loss(u,ww).sqrt())
        assert abs(er-exact)<1e-7,(er,exact)
        if exact<best:best=exact;save(out,u,ww,tree,t,err2=exact)
        row=dict(step=step,error_2=exact,best=best,seconds=time.time()-start,mode=a.mode);trajectory.append(row);print(json.dumps(row),flush=True)
        Path(str(out)+'.json').write_text(json.dumps(trajectory,indent=2))
        if step%100==0:save(HERE/f'trajectories/{out.stem}_{step}.npz',u,ww,tree,t,err2=exact)
