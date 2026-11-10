"""Canonical-angle ADMM refits for the actual two-block PSD/rank-C model.

The unconstrained coefficient step is solved exactly, including interactions
between blocks; the spectral step enforces each block's PSD rank. This offers
a basin-changing weight initialization around unchanged exact integer words.
"""
import argparse,json,time
from core import *
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=3000);p.add_argument('--rho',type=float,default=.01)
p.add_argument('--full-psd',action='store_true');p.add_argument('--fixed-rho',action='store_true');p.add_argument('--warm-dual',action='store_true');a=p.parse_args()
z,u,w,tree,words=load(a.init);R,B,C=w.shape;assert R==2
u=reconstruct_u_vectors(words/128,tree).reshape(R,B,76);sp=np.load(SP);ii,jj,fac=sp['ii'],sp['jj'],sp['fac']
native=HERE/'candidates/native_pair_target.npz'
if not native.exists():
    from native_target import build
    build()
G=np.load(native)['G'];P=u[:,:,ii]*u[:,:,jj]*fac
qs=[];ts=[]
for pp in P:q,tt=np.linalg.qr(pp.T);qs.append(q);ts.append(tt)
left,s,right=np.linalg.svd(qs[0].T@qs[1]);E=qs[0]@left;F=qs[1]@right.T
cc=s[:,None]*s[None,:];h0=E.T@G@E;h1=F.T@G@F;den=1-cc*cc
common=abs(den)<1e-11;AL=np.divide(h0-cc*h1,den,out=h0/2,where=~common);BL=np.divide(h1-cc*h0,den,out=h1/2,where=~common)
D=G-E@AL@E.T-F@BL@F.T;outside=np.square(D).sum()
rank=B if a.full_psd else C
orig=[left.T@ts[0]@w[0],right@ts[1]@w[1]]
Z=np.stack([f@f.T for f in orig]);dual=np.zeros_like(Z);H=np.stack([h0,h1]);rho=a.rho
if a.warm_dual:dual=np.stack([h0-Z[0]-cc*Z[1],h1-Z[1]-cc*Z[0]])/rho
best=float('inf');start=time.time();trajectory=[]
def loss(Z):
    da=Z[0]-AL;db=Z[1]-BL
    return outside+(.5*(1+cc)*np.square(da+db)+.5*(1-cc)*np.square(da-db)).sum()
def save_state(step,err):
    ws=[]
    for r in range(2):
        ev,v=np.linalg.eigh(Z[r]);f=v[:,-rank:]*np.sqrt(np.maximum(ev[-rank:],0))
        ww=np.linalg.solve(ts[r],(left if r==0 else right.T)@f)
        if ww.shape[-1]==C:
            aa,_,bb=np.linalg.svd(ww.T@w[r],full_matrices=False);ww=ww@(aa@bb)
        ws.append(ww)
    ws=np.stack(ws);save(a.out,u,ws,tree,words,err2=err)
    if step%100==0:save(HERE/f'trajectories/{Path(a.out).stem}_{step}.npz',u,ws,tree,words,err2=err)
best=float(np.sqrt(loss(Z)));save_state(-1,best)
for step in range(a.steps):
    rhs=H+rho*(Z-dual);d=1+rho
    A=np.stack([(d*rhs[0]-cc*rhs[1])/(d*d-cc*cc),(d*rhs[1]-cc*rhs[0])/(d*d-cc*cc)])
    old=Z.copy()
    for r in range(2):
        mat=A[r]+dual[r];ev,v=np.linalg.eigh((mat+mat.T)/2);f=v[:,-rank:]*np.sqrt(np.maximum(ev[-rank:],0));Z[r]=f@f.T
    dual+=A-Z
    if step%20==0 or step==a.steps-1:
        er=float(np.sqrt(loss(Z)));primal=float(np.linalg.norm(A-Z));dr=rho*float(np.linalg.norm(Z-old))
        if er<best:best=er;save_state(step,er)
        row=dict(step=step,error_2=er,best=best,rho=rho,primal=primal,dual=dr,rank=rank,linear_relaxation=outside**.5,seconds=time.time()-start)
        trajectory.append(row);print(json.dumps(row),flush=True);Path(str(a.out)+'.json').write_text(json.dumps(trajectory,indent=2))
        if step%100==0 and step>0 and not a.fixed_rho:
            newrho=rho
            if primal>10*dr:newrho=min(1.,rho*2)
            elif dr>10*primal:newrho=max(1e-5,rho/2)
            dual*=rho/newrho;rho=newrho
