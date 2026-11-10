"""Discrete dictionary exchanges using exact projection-span gain screening.

Build new eight-bit vectors from leading residual pair modes; exchange a weak
dictionary row rather than adding rows. The screen computes the exact change
in full dictionary projection energy using its common (B-1)-dimensional span.
Final acceptance refits PSD rank-C weights and checks the actual tensor error.
"""
import argparse,json,time
from core import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=50);p.add_argument('--trials',type=int,default=12);p.add_argument('--seed',type=int,default=5)
a=p.parse_args();torch.manual_seed(a.seed);z,u,w,tree,words=load(a.init);R,B,C=w.shape
words=words.reshape(R,B,75)
m=Metric(tree);ut=m.tensor(u);wt=m.tensor(w);start=time.time()
with torch.no_grad():
 for step in range(a.steps):
    wt,err=m.project(ut,wt,2);r=step%R;Ps=m.feats(ut);P=Ps[r];S=P@P.T;iv=torch.linalg.inv(S)
    Xs=torch.bmm(wt.transpose(1,2),Ps);X=Xs.flatten(0,1);D=m.G-X.T@X
    target=m.G-(X.T@X-Xs[r].T@Xs[r])
    # Randomized residual eigenspace, then exact small Ritz eigensolve.
    Omega=torch.randn((len(D),32),**m.d);Q,_=torch.linalg.qr(D@Omega)
    for _ in range(3):Q,_=torch.linalg.qr(D@Q)
    ev,V=torch.linalg.eigh(Q.T@D@Q);modes=(Q@V[:,-12:]).T
    matrices=torch.zeros((len(modes),76,76),**m.d);matrices[:,m.ii,m.jj]=modes/m.fac
    matrices=matrices+matrices.transpose(1,2)-torch.diag_embed(matrices.diagonal(dim1=-2,dim2=-1))
    ee,vv=torch.linalg.eigh(matrices);ix=ee.abs().argsort(-1,descending=True)[:,:6]
    v=vv.gather(2,ix[:,None,:].expand(-1,76,-1)).transpose(1,2).reshape(-1,76)
    cw=encode(v.cpu().numpy(),tree,8,True);cu=m.decode(m.tensor(cw));F=m.feats(cu)
    # Full-rank projection-gain screen for all B x candidate exchanges.
    dual=iv@P;E=dual/iv.diag().sqrt()[:,None];q=E@F.T
    perp=F-(F@P.T@iv)@P;s=perp.square().sum(-1)
    Z=P@target;ap=Z@perp.T;be=Z@E.T;gp=perp@target
    gpp=(gp*perp).sum(-1);gpe=E@gp.T;gee=(E*(E@target)).sum(-1)
    iva=iv@ap;ivb=iv@be;aa=(ap*iva).sum(0);bb=(be*ivb).sum(0);ab=be.T@iva
    den=s[None,:]+q.square();dproj=(aa[None,:]+2*q*ab+q.square()*bb[:,None]-(gpe+q*gee[:,None]).square())/den.clamp(min=1e-15)
    diag=(gpp[None,:]+2*q*gpe+q.square()*gee[:,None])/den.clamp(min=1e-15)
    gain=2*dproj+diag.square()-(2*bb-gee.square())[:,None]
    gain[den<1e-12]=-torch.inf
    proposals=gain.flatten().argsort(descending=True)[:a.trials].cpu().tolist();changed=False
    for prop in proposals:
        b,j=divmod(prop,len(F));nu=ut.clone();nu[r,b]=cu[j];nw,ne=m.project(nu,wt,2)
        if ne<err-1e-9:
            ut=nu;wt=nw;words[r,b]=cw[j];changed=True;break
    save(a.out,ut,wt,tree,words.reshape(-1,75),err2=ne if changed else err)
    print(json.dumps(dict(step=step,error=err,next=ne,changed=changed,screen_gain=float(gain.max()),seconds=time.time()-start)),flush=True)
    if not changed:break
