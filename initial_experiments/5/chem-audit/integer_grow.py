"""Integer-native residual-directed dictionary growth from selected words.

Every added atom is rounded *before* fitting. Exact block weight refits are
blended with the feasible warm model to retain norm/Lambda limits. This
tests a finite growth proposal family, not an expressivity impossibility.
"""
import argparse,time,json
from audit_core import *
p=argparse.ArgumentParser();p.add_argument('--start',type=int,default=511);p.add_argument('--stop',type=int,default=639)
p.add_argument('--batch',type=int,default=32);p.add_argument('--tag',default='integer_growth');p.add_argument('--lambda-cap',type=float,default=320.);a=p.parse_args()
base=dict(np.load(INC/'submission.npz'));tree=base['tree_edges'];m=Metric(tree);began=time.time()
words=base['angle_words'][76:];u=m.decode(m.tensor(words).reshape(1,719,75));w=m.tensor(m.inp.w);off=m.tensor(m.inp.W)
with torch.no_grad():
    P=m.feats(u)[0];S=P@P.T;K=w[0]@w[0].T;si=torch.linalg.inv(S);diag=si.diag();Y=K@S@K
    importance=2*Y.diag()/diag-(K.diag()/diag).square();ix=importance.argsort(descending=True)[:a.start].sort().values
    words=words[ix.cpu().numpy()];u=u[:,ix];w=w[:,ix]
    # Retain selected physical coefficients initially; a compressed model is
    # not assumed to pass. Its Lambda can be projected toward the cap below.
    rows=[]
    while True:
        B=u.shape[1];cp,_=checkpoint(m,u,w,off,tree,1)
        rows.append(dict(B=B,metrics=cp,seconds=time.time()-began));print(json.dumps(rows[-1]),flush=True)
        save(HERE/f'parameters/{a.tag}_B{B}.npz',u,w,tree,words,W=off.cpu().numpy())
        write_json(HERE/f'results/{a.tag}.json',dict(arguments=vars(a),rows=rows,finished=B>=a.stop))
        if B>=a.stop:break
        P=m.feats(u)[0];F=w[0].T@P;D=m.G-F.T@F
        ev,V=torch.linalg.eigh((D+D.T)/2)
        # Positive residual covariance modes -> full orbital eigenvectors.
        mode=V[:,-12:].T;mat=m.unpack(mode);ee,vv=torch.linalg.eigh(mat)
        cand=vv.transpose(1,2).reshape(-1,76)
        aw=encode(cand.cpu().numpy(),tree,8,True)
        aw=np.unique(aw,axis=0)
        old={tuple(row) for row in words};aw=np.asarray([row for row in aw if tuple(row) not in old])
        uc=m.decode(m.tensor(aw));Pc=m.feats(uc)
        # Score exact decoded atoms outside the current dictionary span.
        Q,_=torch.linalg.qr(P.T);perp=Pc-(Pc@Q)@Q.T
        score=(perp*(perp@D)).sum(-1)/perp.square().sum(-1).clamp(min=1e-16)
        order=score.argsort(descending=True);chosen=[];basis=[]
        for j in order.tolist():
            v=perp[j].clone()
            if basis:
                q=torch.stack(basis);v-=(q@v)@q
            if float(v.norm())<1e-5:continue
            chosen.append(j);basis.append(v/v.norm())
            if len(chosen)>=min(a.batch,a.stop-B):break
        if not chosen:break
        words=np.concatenate([words,aw[chosen]])
        un=torch.cat([u,uc[chosen][None]],dim=1)
        warm=torch.nn.functional.pad(w,(0,0,0,len(chosen)))
        wp,err=m.project(un,warm,1)
        # Largest interpolation step satisfying Lambda/norm caps. Exact
        # residual chooses between the capped refit and zero-weight growth.
        lo,hi=0.,1.
        for _ in range(35):
            alpha=(lo+hi)/2;ww=warm+alpha*(wp-warm)
            g,l,s=m.physical(un,ww,off)
            if float(l)<=a.lambda_cap and float(ww.norm())<=5:lo=alpha
            else:hi=alpha
        trial=warm+lo*(wp-warm)
        w=trial if m.loss(un,trial)<m.loss(un,warm) else warm
        u=un;oldoff=off.clone();m.project_gap(u,w,off)
        if float(m.physical(u,w,off)[1])>a.lambda_cap:off.copy_(oldoff)
