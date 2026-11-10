"""Matched integer-native residual growth versus incumbent-atom reinsertion.

All candidates are quantized before insertion. A common pool and common
exact QR PSD rank-C solves supply relaxed accuracy endpoints; a second
physical-generator transfer records the cost of preserving normalization.
This is an explicit finite proposal family, not a capacity impossibility test.
"""
import argparse
from census_common import *
p=argparse.ArgumentParser();p.add_argument('--start',type=int,default=479);p.add_argument('--stop',type=int,default=511)
p.add_argument('--batch',type=int,default=8);p.add_argument('--tag',default='growth');a=p.parse_args()
z,m,fullu,fullw,off,fullwords=initial();began=time.time()
torch.manual_seed(2026102405)
with torch.no_grad():
    fullP=m.feats(fullu)[0];fullX=fullw[0].T@fullP
    # Use the measured sequential deletion path when available. Otherwise use
    # one-shot leverage for the explicitly separate smaller seed.
    saved=HERE/'parameters'/f'prune_B{a.start}.npz'
    if saved.exists():
        zz,uu,ww,_,aw=load(saved);u=m.tensor(uu);w=m.tensor(ww);words=aw;indices=np.array(zz['indices'])
    else:
        S=fullP@fullP.T;K=fullw[0]@fullw[0].T;iv=torch.linalg.inv(S);d=iv.diag()
        score=2*(K@S@K).diag()/d-(K.diag()/d).square()
        ix=score.argsort(descending=True)[:a.start].sort().values
        indices=ix.cpu().numpy();u=fullu[:,ix];words=fullwords[indices];w,_=m.project(u,fullw[:,ix],1)
    seedu=u.clone();seedw=w.clone();seedwords=words.copy()
    P=m.feats(u)[0];F=w[0].T@P;D=m.G-F.T@F
    # Accurate positive residual Ritz modes; avoid full 2926x2926 eigensolve.
    omega=torch.randn((len(D),48),**m.d);Q,_=torch.linalg.qr(D@omega)
    for _ in range(4):Q,_=torch.linalg.qr(D@Q)
    ev,V=torch.linalg.eigh(Q.T@D@Q);modes=(Q@V[:,-16:]).T
    matrices=m.unpack(modes);ee,vv=torch.linalg.eigh(matrices)
    ix=ee.abs().argsort(-1,descending=True)[:,:8]
    vectors=vv.gather(2,ix[:,None,:].expand(-1,76,-1)).transpose(1,2).reshape(-1,76)
    poolwords=np.unique(encode(vectors.cpu().numpy(),m.tree,8,True),axis=0)
    poolu=m.decode(m.tensor(poolwords));poolP=m.feats(poolu)
    np.savez(HERE/'parameters'/f'{a.tag}_pool.npz',angle_words=poolwords,tree_edges=m.tree)
    rows=[]
    for mode in ['residual','reinsertion']:
        u=seedu.clone();w=seedw.clone();words=seedwords.copy();used=set()
        if mode=='residual':pu=poolu;pw=poolwords;PP=poolP
        else:
            missing=np.setdiff1d(np.arange(719),indices);pu=fullu[0,missing];pw=fullwords[missing];PP=fullP[missing]
        while True:
            B=u.shape[1];P=m.feats(u)[0];Q,_=torch.linalg.qr(P.T)
            transfer=torch.linalg.lstsq(P.T,fullX.T).solution[None]
            oo=off.clone();m.project_gap(u,w,oo)
            ot=off.clone();m.project_gap(u,transfer,ot)
            span=float((m.G-Q@(Q.T@m.G@Q)@Q.T).norm())
            row=dict(mode=mode,B=B,span_bound=span,refit=quantities(m,u,w,oo),transfer=quantities(m,u,transfer,ot),seconds=time.time()-began)
            rows.append(row);print(json.dumps(row),flush=True)
            store(f'{a.tag}_{mode}_B{B}',m,u,w,oo,words)
            store(f'{a.tag}_{mode}_transfer_B{B}',m,u,transfer,ot,words)
            write(f'{a.tag}.json',dict(arguments=vars(a),pool_size=len(pw),rows=rows,complete=mode=='reinsertion' and B>=a.stop))
            if B>=a.stop:break
            F=w[0].T@P;D=m.G-F.T@F
            perp=PP-(PP@Q)@Q.T
            # Score the new direction's exact capture of residual covariance.
            norm=perp.square().sum(-1);v=perp/norm.sqrt().clamp(min=1e-12)[:,None]
            dv=v@D;cross=dv@Q
            gain=2*cross.square().sum(-1)+(dv*v).sum(-1).clamp(min=0).square()
            gain[norm<1e-12]=-torch.inf
            if used:gain[list(used)]=-torch.inf
            order=gain.argsort(descending=True);chosen=[];basis=[]
            for j in order.cpu().tolist():
                if j in used:continue
                vj=perp[j].clone()
                if basis:
                    qb=torch.stack(basis);vj-=(qb@vj)@qb
                if float(vj.norm())<1e-6:continue
                chosen.append(j);basis.append(vj/vj.norm())
                if len(chosen)==min(a.batch,a.stop-B):break
            if not chosen:raise RuntimeError('No independent growth directions')
            used.update(chosen);words=np.concatenate([words,pw[chosen]])
            u=torch.cat([u,pu[chosen][None]],1)
            warm=torch.nn.functional.pad(w,(0,0,0,len(chosen)));w,_=m.project(u,warm,1)
