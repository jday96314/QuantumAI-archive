# Reused from chem-qfit/discrete.py; exact fresh-native Metric via audit_core.
"""Coupled integer-word descent with batched nonlinear row proposals.

Unlike STE or smooth lattice penalties, every evaluated model has exact integer
words. Candidate rows change 1/2/4/8 angles together; their *decoded feature*
differences are screened against the current tensor residual and dictionary
refit tangent space. Independent row proposals are combined, then accepted only
after exact rank-C block weight refits reduce the full signed tensor error.
"""
import argparse,json,time
from audit_core import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=300);p.add_argument('--max-rows',type=int,default=128)
p.add_argument('--cycles',type=int,default=2);p.add_argument('--seed',type=int,default=0)
p.add_argument('--lambda-cap',type=float,default=320.);p.add_argument('--weight-cap',type=float,default=5.)
p.add_argument('--unbounded',action='store_true')
a=p.parse_args();torch.manual_seed(a.seed);z,u,w,tree,words=load(a.init);R,B,C=w.shape
m=Metric(tree);words=encode(u.reshape(-1,76),tree,8,True) if words is None else words
t=m.tensor(np.where(words>=128,words-256,words).reshape(R,B,75));wt=m.tensor(w)
off=m.tensor(z['W']);rejected_by_caps=0;history=[]
def within_caps(uu,ww):
    return a.unbounded or (float(ww.norm())<=a.weight_cap and float(m.physical(uu,ww,off)[1])<=a.lambda_cap)
start=time.time();accepted=0;stale=0
for step in range(a.steps):
    u=m.decode(t);wp,pe=m.project(u,wt,a.cycles)
    if within_caps(u,wp):wt=wp
    else:rejected_by_caps+=1
    err=float(m.loss(u,wt).sqrt())
    save(a.out,u,wt,tree,t,err2=err,W=off.cpu().numpy())
    t.requires_grad_(True);l=m.loss(m.decode(t),wt);l.backward();g=t.grad.detach();t=t.detach()
    with torch.no_grad():
        P=m.feats(u);X=torch.bmm(wt.transpose(1,2),P).flatten(0,1);D=X.T@X-m.G
        Y=torch.bmm(torch.bmm(wt,wt.transpose(1,2)),P);yn=Y.square().sum(-1)
        # Derivative norms identify relevant angles; exact decoded proposals
        # below account for nonlinear rotations and angle coupling.
        tf=t.flatten(0,1);v=torch.zeros((R*B,76),**m.d);v[:,0]=1;J=torch.zeros((R*B,75,76),**m.d)
        for j,(pp,cc) in enumerate(tree):
            co=(tf[:,j]*np.pi/128).cos();sn=(tf[:,j]*np.pi/128).sin()
            jp,jc=J[:,:,pp].clone(),J[:,:,cc].clone();J[:,:,pp]=co[:,None]*jp-sn[:,None]*jc;J[:,:,cc]=sn[:,None]*jp+co[:,None]*jc
            vp,vc=v[:,pp].clone(),v[:,cc].clone();v[:,pp]=co*vp-sn*vc;v[:,cc]=sn*vp+co*vc
            J[:,j,pp]=-np.pi/128*v[:,cc];J[:,j,cc]=np.pi/128*v[:,pp]
        h=(8*J.square().sum(-1).reshape(R,B,75)*yn[...,None]).clamp(min=1e-18)
        # Score each angle by its projected local geometry, not bare gradient.
        for r in range(R):
            S=P[r]@P[r].T;L=torch.linalg.cholesky(S)
            cross=2*(J.reshape(R,B,75,76)[r]@u[r].T)*(u[r]@u[r].T)[:,None,:]
            cr=cross.reshape(B*75,B);proj=(cr*torch.cholesky_solve(cr.T,L).T).sum(-1).reshape(B,75)
            raw=2*J.reshape(R,B,75,76)[r].square().sum(-1)
            h[r]=4*((raw-proj).clamp(min=0)+.001*raw)*yn[r,:,None]
        scores=.5*h-g.abs();scores[h<1e-16]=torch.inf
        order=scores.argsort(-1);direction=-g.sign();candidates=[];ident=[];changes=[]
        for r in range(R):
            cp=[]
            for k in [1,2,4,8]:
                ix=order[r,:,:k];delta=torch.zeros_like(t[r]);delta.scatter_(1,ix,direction[r].gather(1,ix))
                cp.append(t[r]+delta)
            cp=torch.stack(cp,1);dp=m.feats(m.decode(cp))-P[r,:,None,:]
            flat=dp.reshape(-1,dp.shape[-1]);cross=flat@P[r].T
            S=P[r]@P[r].T;L=torch.linalg.cholesky(S)
            perp=(flat.square().sum(-1)-(cross*torch.cholesky_solve(cross.T,L).T).sum(-1)).clamp(min=0)
            linear=4*(dp*(Y[r]@D)[:,None,:]).sum(-1)
            pred=linear+2*yn[r,:,None]*perp.reshape(B,4)
            pred[dp.square().sum(-1)<1e-22]=torch.inf
            best,ki=pred.min(-1)
            for b in range(B):
                candidates.append(float(best[b]));ident.append((r,b));changes.append(cp[b,ki[b]]-t[r,b])
        ix=np.argsort(candidates);groups=[]
        # Include coupled row moves even when no individual fixed-weight move
        # passes; the exact weight refit resolves the interaction.
        for k in [a.max_rows,64,32,16,8,4,2,1]:
            k=min(k,R*B)
            if k>0:groups.append(ix[:k])
        groups.extend([ix[j:j+1] for j in range(1,12)])
        changed=False
        for group in groups:
            trial=t.clone()
            for j in group:
                r,b=ident[j];trial[r,b]+=changes[j]
            uu=m.decode(trial);ww,ne=m.project(uu,wt,a.cycles)
            feasible=within_caps(uu,ww)
            if not feasible:rejected_by_caps+=1
            if ne<err-1e-10 and feasible:
                t=trial;wt=ww;accepted+=len(group);changed=True;break
        stale=0 if changed else stale+1
    if step%5==0 or not changed:
        cp,_=checkpoint(m,m.decode(t),wt,off,tree,1)
        row=dict(step=step,error=err,next=ne,accepted_rows=accepted,changed=changed,metrics=cp,rejected_by_caps=rejected_by_caps,s=time.time()-start)
        history.append(row);print(json.dumps(row),flush=True)
    if stale>=3:break
with torch.no_grad():
    u=m.decode(t);wp,pe=m.project(u,wt,6)
    if within_caps(u,wp):wt=wp
    err=float(m.loss(u,wt).sqrt());oldoff=off.clone();m.project_gap(u,wt,off)
    if not within_caps(u,wt):off.copy_(oldoff)
    save(a.out,u,wt,tree,t,err2=err,W=off.cpu().numpy())
    cp,_=checkpoint(m,u,wt,off,tree,2)
rec=dict(arguments=vars(a),final=err,metrics=cp,history=history,rejected_by_caps=rejected_by_caps,seconds=time.time()-start)
write_json(HERE/f'results/{Path(a.out).stem}.json',rec);print(json.dumps(rec),flush=True)
