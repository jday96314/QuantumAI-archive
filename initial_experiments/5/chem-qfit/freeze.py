"""Progressive vector freezing with tensor-error feedback and residual retraining.

Train remaining continuous word coordinates and all weights against the actual
tensor residual after every frozen batch. Frozen vectors stay exactly integer;
the final saved model is decoded and refitted in float64, with exact words.
"""
import argparse,json,time
from core import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--warmup',type=int,default=500);p.add_argument('--steps',type=int,default=150)
p.add_argument('--batch',type=int,default=40);p.add_argument('--lr',type=float,default=.06)
p.add_argument('--weight-lr',type=float,default=1e-4);p.add_argument('--seed',type=int,default=0)
p.add_argument('--order',choices=['feedback','random','hardest'],default='feedback');p.add_argument('--double',action='store_true')
p.add_argument('--unit',choices=['vector','angle'],default='vector')
p.add_argument('--probe-only',action='store_true')
p.add_argument('--no-align',action='store_true',help='matched gauge-continuity ablation')
a=p.parse_args();torch.manual_seed(a.seed);z,u,w,tree,words=load(a.init);R,B,C=w.shape
m=Metric(tree,torch.float64 if a.double else torch.float32);md=Metric(tree)
md.align_columns=not a.no_align
if words is None:
    v=u.reshape(-1,76).copy();th=np.zeros((R*B,75))
    for j in reversed(range(75)):
        pp,cc=tree[j];sg=np.where(v[:,pp]<0,-1.,1.)
        th[:,j]=np.arctan2(v[:,cc]*sg,v[:,pp]*sg)*128/np.pi
        v[:,pp]=sg*np.hypot(v[:,pp],v[:,cc]);v[:,cc]=0
else:th=np.where(words>=128,words-256,words).astype(float)
t=m.tensor(th.reshape(R,B,75)).requires_grad_(True);wt=m.tensor(w).requires_grad_(True)
frozen=torch.zeros(t.shape if a.unit=='angle' else (R,B),device='cuda',dtype=torch.bool);fixed=t.detach().round().clone();start=time.time()
opt=torch.optim.Adam([{'params':[t],'lr':a.lr},{'params':[wt],'lr':a.weight_lr}])
def training(steps):
    for step in range(steps):
        opt.zero_grad();l=m.loss(m.decode(t),wt);l.backward()
        t.grad[frozen]=0;opt.step()
        with torch.no_grad():t[frozen]=fixed[frozen]
    return float(m.loss(m.decode(t),wt).detach().sqrt())
print(json.dumps(dict(phase='warmup',error=training(a.warmup),s=time.time()-start)),flush=True)
with torch.no_grad():
    ud=md.decode(t.detach().double());wd,soft=md.project(ud,wt.detach().double(),2)
    uq=md.decode(t.detach().double().round());direct=float(md.loss(uq,wd).sqrt());wq,hard=md.project(uq,wd,4)
    prefix=str(Path(a.out).with_suffix(''))
    save(prefix+'_warmup_float.npz',ud,wd,tree,err2=soft)
    save(prefix+'_warmup_integer.npz',uq,wq,tree,t.round(),err2=hard)
    print(json.dumps(dict(phase='matched_probe',continuous=soft,direct=direct,refitted=hard,s=time.time()-start)),flush=True)
if a.probe_only:raise SystemExit(0)
best=1e9
for phase in range((frozen.numel()+a.batch-1)//a.batch):
    with torch.no_grad():
        uu=m.decode(t);uq=m.decode(t.round());P=m.feats(uu);dp=m.feats(uq)-P
        X=torch.bmm(wt.transpose(1,2),P).flatten(0,1);D=X.T@X-m.G
        # Leading exact residual change plus the norm of the rank-two local
        # tensor perturbation. This ranks *tensor* impact, not angle distance.
        Y=torch.bmm(torch.bmm(wt,wt.transpose(1,2)),P);yn=Y.square().sum(-1)
        dd=dp.square().sum(-1);yd=(Y*dp).sum(-1);wa=wt.square().sum(-1)
        score=4*(dp*(Y@D)).sum(-1)+2*wa*((dp@D)*dp).sum(-1)+2*yn*dd+2*yd.square()+4*wa*yd*dd+wa.square()*dd.square()
        if a.unit=='angle':
            # Tensor Gauss-Newton coordinate curvature. Small-angle freezes permit compensation
            # by the other angles in the same vector.
            tf=t.flatten(0,1);v=torch.zeros((R*B,76),**m.d);v[:,0]=1
            J=torch.zeros((R*B,75,76),**m.d)
            for j,(pp,cc) in enumerate(tree):
                co=(tf[:,j]*np.pi/128).cos();sn=(tf[:,j]*np.pi/128).sin()
                jp,jc=J[:,:,pp].clone(),J[:,:,cc].clone();J[:,:,pp]=co[:,None]*jp-sn[:,None]*jc;J[:,:,cc]=sn[:,None]*jp+co[:,None]*jc
                vp,vc=v[:,pp].clone(),v[:,cc].clone();v[:,pp]=co*vp-sn*vc;v[:,cc]=sn*vp+co*vc
                J[:,j,pp]=-np.pi/128*v[:,cc];J[:,j,cc]=np.pi/128*v[:,pp]
            raw=2*J.square().sum(-1).reshape(R,B,75)
            h=4*raw*yn[...,None]
            delta=t.round()-t
        if a.order=='random':score=torch.rand_like(t if a.unit=='angle' else score)
    if a.unit=='angle' and a.order!='random':
        opt.zero_grad();l=m.loss(m.decode(t),wt);l.backward()
        with torch.no_grad():score=t.grad*delta+.5*h*delta.square()
    with torch.no_grad():
        if a.order=='hardest':score=-score
        score[frozen]=torch.inf;count=min(a.batch,int((~frozen).sum()));ix=score.flatten().argsort()[:count]
        frozen.flatten()[ix]=True
        if a.unit=='angle':fixed.flatten()[ix]=t.flatten()[ix].round()
        else:fixed.flatten(0,1)[ix]=t.flatten(0,1)[ix].round()
        t[frozen]=fixed[frozen]
    err=training(a.steps)
    with torch.no_grad():
        th=t.detach().double();ud=md.decode(th);wd,continuous=md.project(ud,wt.detach().double(),2)
        wt.copy_(wd.to(wt.dtype));uq=md.decode(th.round());wq,hard=md.project(uq,wd,4)
        if hard<best:
            best=hard;save(a.out,uq,wq,tree,th.round(),err2=hard)
        save(str(Path(a.out).with_suffix(''))+'_soft.npz',ud,wd,tree,err2=continuous,theta=th.cpu().numpy(),frozen=frozen.cpu().numpy())
    print(json.dumps(dict(phase=phase,frozen=int(frozen.sum()),soft=continuous,hard=hard,best=best,s=time.time()-start)),flush=True)
