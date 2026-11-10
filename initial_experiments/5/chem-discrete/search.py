"""Cross-block finite-feature QUBO and true integer-temperature basin search.

Every state and proposal is an integer word table. Residual Ritz modes generate
correlated moves; finite feature differences and cross-block Hessian terms
screen bundles. Acceptance uses the full pair tensor with separate PSD rank-C
block refits. No continuous vector fitting or angle rounding is performed.
"""
import argparse,json,time
from compact import *
from lattice_qubo import solve as qubo

p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=300);p.add_argument('--seconds',type=float,default=1200)
p.add_argument('--seed',type=int,default=1);p.add_argument('--pool',type=int,default=384)
p.add_argument('--mode',choices=['coupled','thermal','single'],default='coupled')
p.add_argument('--cycles',type=int,default=3);p.add_argument('--temperature',type=float,default=1e-7)
p.add_argument('--refit-budget',type=int,default=0);p.add_argument('--norm-weight',type=float,default=0.)
p.add_argument('--subtree',action='store_true')
p.add_argument('--screen32',action='store_true')
p.add_argument('--legacy-temperature',action='store_true')
p.add_argument('--audit-rank-init',default='')
a=p.parse_args();torch.manual_seed(a.seed);rng=np.random.default_rng(a.seed)
z,u,w,tree,words=load(a.init);R,B,C=w.shape;m=Metric(tree)
t=m.tensor(words.reshape(R,B,75));wt=m.tensor(w);u=m.decode(t)
wt,err=m.project(u,wt,20);best=err;bestt=t.clone();bestw=wt.clone();start=time.time();out=Path(a.out)
save(out,u,wt,tree,t,err2=err);trajectory=[];stale=0;uphill=0;accepted=0
refits=20*R;norm=None;normbest=float('inf');off=None;normrec={}
if a.norm_weight:
    from normalization import Normalization
    off=z.get('W',np.zeros((R,C)));norm=Normalization(m,off)
    normvalue,off,normrec=norm.evaluate(u,wt,40);norm.off=off
    merit=err*err+a.norm_weight*normvalue
bestoff=off.clone() if off is not None else None
bestmerit=merit if norm else None;bestnormrec=normrec.copy()
auditw=m.tensor(load(a.audit_rank_init)[2]) if a.audit_rank_init else None
auditbest=float('inf');auditrecords=[]

@torch.no_grad()
def propose(t,w,step):
    u,J=derivatives(m,t);P=m.feats(u);X=torch.bmm(w.transpose(1,2),P).flatten(0,1)
    D=X.T@X-m.G;Y=torch.bmm(torch.bmm(w,w.transpose(1,2)),P)
    # Use several signed physical residual modes, rather than a span-only loss.
    dused=D
    if step%4 in [1,2]:
        omega=torch.randn((2926,12),**m.d);q=torch.linalg.qr(D@omega).Q
        for _ in range(2):q=torch.linalg.qr(D@q).Q
        ev,v=torch.linalg.eigh(q.T@D@q);ids=ev.abs().argsort(descending=True)[:(2 if step%4==1 else 6)]
        modes=q@v[:,ids];dused=(modes*ev[ids])@modes.T
    yd=Y@dused
    matrices=torch.zeros((R,B,76,76),**m.d);matrices[:,:,m.ii,m.jj]=yd*m.fac/2
    matrices=matrices+matrices.transpose(-1,-2)
    # Diagonal was doubled in the preceding symmetric construction.
    g=8*(J*(matrices@u[...,None]).squeeze(-1)[:,:,None,:]).sum(-1)
    raw=8*J.square().sum(-1)*Y.square().sum(-1)[...,None]
    # Projection geometry for the most relevant candidate derivatives only.
    rankscore=g.abs()/raw.sqrt().clamp(min=1e-12)
    order=rankscore.argsort(-1,descending=True)
    cwords=[];rows=[]
    for k in ([1] if a.mode=='single' else [1,2,4,8]):
        inds=order[:,:,:k];delta=torch.zeros_like(t)
        delta.scatter_(-1,inds,-g.gather(-1,inds).sign())
        if step%7==6 and a.mode!='single':delta*=2
        cwords.append((t+delta).reshape(R*B,75));rows.append(torch.arange(R*B,device='cuda'))
    # Random local tree transformations supply directions cancelled by g-sum.
    if a.mode!='single':
        inds=order[:,:,:8];scores=torch.rand_like(inds,dtype=torch.float64)
        use=scores>.55;delta=torch.zeros_like(t)
        ss=torch.where(torch.rand_like(scores)>.2,-g.gather(-1,inds).sign(),g.gather(-1,inds).sign())
        delta.scatter_(-1,inds,ss*use);cwords.append((t+delta).reshape(R*B,75));rows.append(torch.arange(R*B,device='cuda'))
    if a.subtree:
        # Mutate complete short Givens branches. Internal zero-support
        # angles have no local derivative; a root activation plus a branch
        # replacement can reach directions unavailable to coordinate descent.
        desc={j:{int(c)} for j,(p,c) in enumerate(tree)}
        for j in range(74,-1,-1):
            for k,(p,c) in enumerate(tree):
                if k>j and int(p) in desc[j]:desc[j].add(int(c))
        branches=[(j,[k for k,(p,c) in enumerate(tree) if k>j and int(p) in ds]) for j,ds in desc.items() if 2<=len(ds)<=20]
        for repeat in range(3):
            changed=t.clone().reshape(R*B,75)
            for j,ks in branches:
                chosen=torch.rand((R*B,),device='cuda')<1/max(1,len(branches))
                rr0=torch.nonzero(chosen,as_tuple=True)[0]
                if len(rr0):
                    changed[rr0[:,None],torch.tensor(ks,device='cuda')[None,:]]+=torch.randint(-32,33,(len(rr0),len(ks)),device='cuda')
                    changed[rr0,j]+=torch.where(torch.rand((len(rr0),),device='cuda')>.5,1.,-1.)
            cwords.append(changed);rows.append(torch.arange(R*B,device='cuda'))
    cw=torch.cat(cwords);rr=torch.cat(rows);pf=P.flatten(0,1);yf=Y.flatten(0,1)
    dp=m.feats(m.decode(cw))-pf[rr];perp=dp.float() if a.screen32 else dp.clone()
    for r in range(R):
        ix=(rr//B==r);q=torch.linalg.qr(P[r].T).Q
        q=q.to(perp.dtype);perp[ix]-=(perp[ix]@q)@q.T
    ys=yf[rr].to(perp.dtype);lin=4*(perp*(ys@D.to(perp.dtype))).sum(-1)
    diag=4*(perp.square().sum(-1)*ys.square().sum(-1)+(perp*ys).sum(-1).square())
    unary=lin+.5*diag
    unary[dp.square().sum(-1)<1e-20]=torch.inf
    # Keep multiple proposal families represented, including individually
    # uphill candidates that may help a coordinated donor/receiver move.
    ids=unary.argsort()[:a.pool]
    if a.mode!='single':
        ids=torch.unique(torch.cat([ids,(unary/diag.clamp(min=1e-16)).argsort()[:a.pool//3]]))[:a.pool]
    cw,rr,perp,ys,lin,diag,unary=[x[ids] for x in [cw,rr,perp,ys,lin,diag,unary]]
    h=4*((perp@perp.T)*(ys@ys.T)+(perp@ys.T)*(ys@perp.T))
    # Same-row alternative candidates are mutually exclusive, not additive.
    conflict=rr[:,None]==rr[None,:];h[conflict]=torch.maximum(h[conflict],torch.tensor(1e3,**m.d));h.fill_diagonal_(0)
    proposals=[]
    if a.mode!='single':
        for cap in [16,64,192]:
            chosen,energy=qubo(h.cpu().numpy(),unary.cpu().numpy(),restarts=8,cap=cap,seed=a.seed*100000+step*11+cap)
            if len(chosen):proposals.append((chosen,energy))
    order=unary.argsort().cpu().numpy();seen=set();distinct=[]
    for j in order:
        row=int(rr[j])
        if row not in seen:distinct.append(j);seen.add(row)
    for count in ([1] if a.mode=='single' else [64,16,4,1]):
        js=np.array(distinct[:count]);proposals.append((js,float(unary[js].sum())))
    if a.mode=='thermal':
        for _ in range(3):
            js=rng.choice(len(cw),size=min(4,len(cw)),replace=False)
            if len(set(rr[js].cpu().tolist()))==len(js):proposals.append((js,float(unary[js].sum())))
    return cw,rr,proposals

for step in range(a.steps):
    cw,rr,proposals=propose(t,wt,step);changed=False;ne=err;family=None
    if a.mode=='thermal' and not a.legacy_temperature and step%7==6:rng.shuffle(proposals)
    for pi,(chosen,pred) in enumerate(proposals):
        trial=t.clone().reshape(R*B,75);trial[rr[chosen]]=cw[chosen];trial=trial.reshape(R,B,75)
        uu=m.decode(trial);ww,ne=m.project(uu,wt,a.cycles);refits+=a.cycles*R
        if norm:
            nv,no,nr=norm.evaluate(uu,ww);newmerit=ne*ne+a.norm_weight*nv
        temp=a.temperature*(.1+.9*(1-(step%50)/50))
        thermal_delta=newmerit-merit if norm else ne*ne-err*err
        thermal=(a.mode=='thermal' and (stale>2 or (not a.legacy_temperature and step%7==6)) and thermal_delta<20*temp and rng.random()<np.exp(min(0,-thermal_delta/temp)))
        if (newmerit<merit-1e-14 if norm else ne<err-1e-10) or thermal:
            uphill+=int(ne>err);t,wt,err=trial,ww,ne;changed=True;accepted+=len(chosen);family=pi;break
        if a.refit_budget and refits>=a.refit_budget:break
    if norm and changed:
        merit=newmerit;off=no;norm.off=off;normrec=nr
        if merit<normbest:
            normbest=merit;save(str(out.with_suffix(''))+'_merit.npz',m.decode(t),wt,tree,t,W=off.cpu().numpy(),err2=err)
    if err<best:
        best=err;bestt=t.clone();bestw=wt.clone()
        if norm:bestoff=off.clone();bestmerit=merit;bestnormrec=normrec.copy()
        kw=dict(W=off.cpu().numpy()) if norm else {}
        save(out,m.decode(t),wt,tree,t,err2=err,**kw)
        save(HERE/f'trajectories/{out.stem}_best{step}.npz',m.decode(t),wt,tree,t,err2=err,**kw)
    stale=0 if changed and err<=best else stale+1
    if a.mode=='thermal' and step%50==49:
        t=bestt.clone();wt=bestw.clone();err=best;stale=0
        if norm:off=bestoff.clone();norm.off=off;merit=bestmerit;normrec=bestnormrec.copy()
    row=dict(step=step,error_2=err,best=best,changed=changed,accepted_words_or_rows=accepted,uphill=uphill,proposal_family=family,PSD_block_sweeps=refits,normalization=normrec,metric_target='native_float32_reference',seconds=time.time()-start)
    trajectory.append(row);print(json.dumps(row),flush=True)
    Path(str(out)+'.json').write_text(json.dumps(trajectory,indent=2))
    if auditw is not None and (step%20==0 or step==a.steps-1 or time.time()-start>=a.seconds):
        auditw,ae=m.project(m.decode(bestt),auditw,20)
        ar=dict(step=step,R=R,B=B,C=auditw.shape[-1],error_2=ae,geometry_error_2=best,seconds=time.time()-start)
        auditrecords.append(ar);print(json.dumps(dict(rank_audit=ar)),flush=True)
        save(HERE/f'trajectories/{out.stem}_rank{auditw.shape[-1]}_{step}.npz',m.decode(bestt),auditw,tree,bestt,err2=ae)
        if ae<auditbest:
            auditbest=ae;save(HERE/f'candidates/{out.stem}_rank{auditw.shape[-1]}_best.npz',m.decode(bestt),auditw,tree,bestt,err2=ae)
        Path(str(out)+'_rank_audit.json').write_text(json.dumps(auditrecords,indent=2))
    if step%20==0:save(HERE/f'trajectories/{out.stem}_current{step}.npz',m.decode(t),wt,tree,t,err2=err)
    if time.time()-start>=a.seconds:break
    if a.refit_budget and refits>=a.refit_budget:break
print(json.dumps(dict(final=best,seconds=time.time()-start,mode=a.mode)),flush=True)
