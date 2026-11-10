"""Validate proposal derivatives against the true signed physical tensor."""
import json
from compact import *
z,u,w,tree,words=load(HERE/'candidates/start_legacy_exact_split_B479.npz');m=Metric(tree);R,B,C=w.shape
t=m.tensor(words.reshape(R,B,75));wt=m.tensor(w);th=t.clone().requires_grad_(True)
loss=m.loss(m.decode(th),wt);auto=torch.autograd.grad(loss,th)[0]
with torch.no_grad():
    u,J=derivatives(m,t);P=m.feats(u);X=torch.bmm(wt.transpose(1,2),P).flatten(0,1);D=X.T@X-m.G;Y=torch.bmm(wt@wt.transpose(1,2),P);yd=Y@D
    M=torch.zeros((R,B,76,76),**m.d);M[:,:,m.ii,m.jj]=yd*m.fac/2;M=M+M.transpose(-1,-2)
    manual=8*(J*(M@u[...,None]).squeeze(-1)[:,:,None,:]).sum(-1)
    checks=[]
    for index in auto.abs().flatten().argsort(descending=True)[:8].cpu().tolist():
        plus=t.clone();minus=t.clone();plus.flatten()[index]+=.001;minus.flatten()[index]-=.001
        fd=float((m.loss(m.decode(plus),wt)-m.loss(m.decode(minus),wt))/.002)
        checks.append(dict(index=index,autograd=float(auto.flatten()[index]),manual=float(manual.flatten()[index]),finite_difference=fd))
    row=dict(max_analytic_difference=float((manual-auto).abs().max()),max_finite_difference_error=max(abs(x['autograd']-x['finite_difference']) for x in checks),checks=checks)
    assert row['max_analytic_difference']<1e-10 and row['max_finite_difference_error']<1e-9,row
(HERE/'results/proposal_gradient_control.json').write_text(json.dumps(row,indent=2));print(json.dumps(row))
