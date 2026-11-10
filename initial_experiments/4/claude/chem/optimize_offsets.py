"""Minimize the convex SOS gap for fixed X, then choose an accepted offset."""
import argparse,json,time
from pathlib import Path
import numpy as np
import torch
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=80);a=p.parse_args();torch.set_num_threads(3)
HERE=Path(__file__).resolve().parent;z=dict(np.load(a.fit));ref=np.load(Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters')/'chem_reference.npz')
u=z['U'];w=z['w'];R,B,C=w.shape;n=76;d=dict(device='cuda',dtype=torch.float64)
X=np.einsum('rbp,rbq,rbc->rcpq',u,u,w,optimize=True).reshape(-1,n,n)
Xt=torch.tensor(X,**d);wt=torch.tensor(w,**d)
h0=torch.tensor(ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(n),**d)
constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
W=torch.zeros(R*C,**d,requires_grad=True)
def vals(off):
    h=h0-torch.einsum('c,cpq->pq',off,Xt);eig=torch.linalg.eigvalsh(.5*(h+h.T))
    norm=torch.sqrt(eig.square()+1e-14).sum()
    gap=constant+norm-eig.sum()+.5*off.square().sum()
    wb=off.reshape(R,C)+wt.sum(1)
    lam=norm+.25*(wt.abs().sum(1)+wb.abs()).square().sum()
    return gap,lam
opt=torch.optim.LBFGS([W],lr=1,max_iter=12,history_size=30,line_search_fn='strong_wolfe',tolerance_grad=1e-9,tolerance_change=1e-13)
start=time.time()
for s in range(a.steps):
    def closure():
        opt.zero_grad();g,_=vals(W);g.backward();return g
    opt.step(closure)
    if s%10==0:
        g,l=vals(W);print(json.dumps(dict(step=s,gap=float(g),Lambda=float(l),seconds=time.time()-start)),flush=True)
with torch.no_grad():
    g,l=vals(W)
    if float(g)<3.916:
        origin=W.clone();direction=-wt.sum(1).reshape(-1)-origin
        lo,hi=0.,1.
        if float(vals(origin+direction)[0])<3.916:direction=-origin
        for _ in range(80):
            mid=(lo+hi)/2
            if float(vals(origin+mid*direction)[0])<3.916:lo=mid
            else:hi=mid
        W.copy_(origin+hi*direction)
    g,l=vals(W);print(json.dumps(dict(final_gap=float(g),Lambda=float(l),lambda_eff=float((g*(2*l-g)).sqrt()))),flush=True)
    z['W']=W.cpu().numpy().reshape(R,C)
np.savez(a.out,**z)
