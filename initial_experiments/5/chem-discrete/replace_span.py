"""Integer atom-bank exchange, coordinated block donors, and thermal escapes.

The relaxed span energy only orders trials. Every accepted state has a full
signed pair-tensor residual and separate exact PSD rank-C block refits.
"""
import argparse,json,time
from core import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=200);p.add_argument('--seconds',type=float,default=1000)
p.add_argument('--seed',type=int,default=13);p.add_argument('--trials',type=int,default=12)
p.add_argument('--temperature',type=float,default=2e-7);a=p.parse_args()
torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed);z,u,w,tree,words=load(a.init)
R,B,C=w.shape;m=Metric(tree);t=m.tensor(words.reshape(R,B,75));ut=m.decode(t);wt=m.tensor(w)
wt,err=m.project(ut,wt,20);best=err;bt=t.clone();bw=wt.clone();start=time.time();out=Path(a.out)
save(out,ut,wt,tree,t,err2=err);trajectory=[];uphill=0;accepted=0
with torch.no_grad():
 for step in range(a.steps):
    r=step%R;Pall=m.feats(ut);P=Pall[r];Xs=torch.bmm(wt.transpose(1,2),Pall);X=Xs.flatten(0,1)
    D=m.G-X.T@X;target=D+Xs[r].T@Xs[r]
    S=P@P.T;iv=torch.linalg.inv(S)
    # Candidate bank contains exact other-block words and entirely new
    # representable atoms suggested by signed physical residual modes.
    bank=[]
    if R>1:bank.append(t[(r+1)%R])
    omega=torch.randn((2926,24),**m.d);q=torch.linalg.qr(D@omega).Q
    for _ in range(3):q=torch.linalg.qr(D@q).Q
    ev,v=torch.linalg.eigh(q.T@D@q);ids=ev.abs().argsort(descending=True)[:10];modes=(q@v[:,ids]).T
    matrices=torch.zeros((len(ids),76,76),**m.d);matrices[:,m.ii,m.jj]=modes/m.fac
    matrices=matrices+matrices.transpose(1,2)-torch.diag_embed(matrices.diagonal(dim1=-2,dim2=-1))
    ee,vv=torch.linalg.eigh(matrices);ix=ee.abs().argsort(-1,descending=True)[:,:8]
    vectors=vv.gather(2,ix[:,None,:].expand(-1,76,-1)).transpose(1,2).reshape(-1,76)
    bank.append(m.tensor(encode(vectors.cpu().numpy(),tree,8,True)))
    # Quantized sums/differences are dictionary proposals only, never float
    # reconstruction fits. Their integer words are evaluated immediately.
    pairs=torch.randint(B,(48,2),device='cuda');vv=ut[r,pairs[:,0]]+torch.where(torch.rand((48,1),**m.d)>.5,1.,-1.)*ut[r,pairs[:,1]]
    vv=vv/vv.norm(dim=-1,keepdim=True).clamp(min=1e-20)
    bank.append(m.tensor(encode(vv.cpu().numpy(),tree,8,True)))
    cw=torch.cat(bank);cu=m.decode(cw);F=m.feats(cu)
    dual=iv@P;E=dual/iv.diag().sqrt()[:,None];qf=E@F.T
    perp=F-(F@P.T@iv)@P;s=perp.square().sum(-1)
    Z=P@target;ap=Z@perp.T;be=Z@E.T;gp=perp@target
    gpp=(gp*perp).sum(-1);gpe=E@gp.T;gee=(E*(E@target)).sum(-1)
    iva=iv@ap;ivb=iv@be;aa=(ap*iva).sum(0);bb=(be*ivb).sum(0);ab=be.T@iva
    den=s[None,:]+qf.square()
    dproj=(aa[None,:]+2*qf*ab+qf.square()*bb[:,None]-(gpe+qf*gee[:,None]).square())/den.clamp(min=1e-15)
    diag=(gpp[None,:]+2*qf*gpe+qf.square()*gee[:,None])/den.clamp(min=1e-15)
    gain=2*dproj+diag.square()-(2*bb-gee.square())[:,None];gain[den<1e-10]=-torch.inf
    choices=gain.flatten().argsort(descending=True)[:max(64,a.trials)].cpu().tolist()
    distinct=[];used=set()
    for prop in choices:
        b,j=divmod(prop,len(F))
        if b not in used:distinct.append((b,j));used.add(b)
    proposals=[distinct[:k] for k in [4,2] if len(distinct)>=k]+[[(p//len(F),p%len(F))] for p in choices[:a.trials]]
    changed=False;ne=err
    for pi,proposal in enumerate(proposals):
        trial=t.clone()
        for b,j in proposal:trial[r,b]=cw[j]
        # In donor/receiver proposals, also move the displaced integer atom
        # into its donor location; this changes the two block spans together.
        if step%3==2 and R>1:
            for b,j in proposal:
                if j<B:trial[(r+1)%R,j]=t[r,b]
        uu=m.decode(trial);ww,ne=m.project(uu,wt,4)
        temp=a.temperature*(.1+.9*(1-(step%40)/40))
        thermal=ne*ne-err*err<15*temp and rng.random()<np.exp(min(0,(err*err-ne*ne)/temp))
        if ne<err-1e-10 or thermal:
            uphill+=int(ne>err);t,ut,wt,err=trial,uu,ww,ne;accepted+=len(proposal);changed=True;break
    if err<best:
        best=err;bt=t.clone();bw=wt.clone();save(out,ut,wt,tree,t,err2=err)
        save(HERE/f'trajectories/{out.stem}_best{step}.npz',ut,wt,tree,t,err2=err)
    if step%40==39:t=bt.clone();wt=bw.clone();ut=m.decode(t);err=best
    row=dict(step=step,error_2=err,best=best,changed=changed,accepted_rows=accepted,uphill=uphill,proposal_family=pi,screen_gain=float(gain.max()),seconds=time.time()-start)
    trajectory.append(row);print(json.dumps(row),flush=True);Path(str(out)+'.json').write_text(json.dumps(trajectory,indent=2))
    if time.time()-start>=a.seconds:break
