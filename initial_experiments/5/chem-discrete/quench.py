"""Short exact integer quench between combinatorial basin changes."""
from compact import derivatives
from core import *

def quench(m,t,w,steps=4):
    R,B,_=t.shape;count=0
    for step in range(steps):
        with torch.enable_grad():
            th=t.clone().requires_grad_(True);loss=m.loss(m.decode(th),w);g=torch.autograd.grad(loss,th)[0]
        with torch.no_grad():
            u,J=derivatives(m,t);P=m.feats(u);Y=torch.bmm(w@w.transpose(1,2),P)
            h=(8*J.square().sum(-1)*Y.square().sum(-1)[...,None]).clamp(min=1e-18)
            relevance=g.abs()/h.sqrt();ix=relevance.argmax(-1)
            delta=torch.zeros_like(t);delta.scatter_(-1,ix[...,None],-g.gather(-1,ix[...,None]).sign())
            cand=t+delta;dp=m.feats(m.decode(cand))-P;X=torch.bmm(w.transpose(1,2),P).flatten(0,1);D=X.T@X-m.G
            perp=dp.clone()
            for r in range(R):
                Q=torch.linalg.qr(P[r].T).Q;perp[r]-=(dp[r]@Q)@Q.T
            score=4*(perp*(Y@D)).sum(-1)+2*Y.square().sum(-1)*perp.square().sum(-1)
            order=score.flatten().argsort();err=float(loss.detach().sqrt());changed=False
            for k in [64,16,4,1]:
                ids=order[:k];trial=t.clone().reshape(R*B,75);trial[ids]=cand.reshape(R*B,75)[ids];trial=trial.reshape_as(t)
                ww,ne=m.project(m.decode(trial),w,2);count+=2
                if ne<err-1e-10:t,w=trial,ww;changed=True;break
            if not changed:break
    return t,w,count
