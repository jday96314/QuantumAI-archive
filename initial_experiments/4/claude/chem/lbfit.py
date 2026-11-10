"""L-BFGS fit of a block DF-THC shape with periodic exact w projection."""
import argparse, json, time
import numpy as np, torch
p=argparse.ArgumentParser()
for k in 'R B C'.split():p.add_argument('--'+k,type=int,required=True)
p.add_argument('--outer',type=int,default=200);p.add_argument('--inner',type=int,default=50)
p.add_argument('--init',required=True);p.add_argument('--out',required=True);p.add_argument('--cycles',type=int,default=2)
a=p.parse_args()
sp=np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev,V=sp['eigenvalues'],sp['eigenvectors']
dt=torch.float64;d=dict(device='cuda',dtype=dt)
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d)
ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
R,B,C,n=a.R,a.B,a.C,76
def feats(u):
    un=u/u.norm(dim=-1,keepdim=True);return un[...,ii]*un[...,jj]*fac
def project(u,w,cycles):
    with torch.no_grad():
        F=feats(u);cont=[]
        for r in range(R):
            Xr=w[r].T@F[r];cont.append(Xr.T@Xr)
        tot=sum(cont)
        for _ in range(cycles):
            for r in range(R):
                Gr=G-(tot-cont[r]);Q,Rm=torch.linalg.qr(F[r].T)
                Hp=Q.T@Gr@Q;Hp=(Hp+Hp.T)/2;e,v=torch.linalg.eigh(Hp)
                idx=torch.argsort(e,descending=True)[:C];root=v[:,idx]*e[idx].clamp(min=0).sqrt()
                wr=torch.linalg.solve(Rm,root);w[r]=wr
                Xr=wr.T@F[r];new=Xr.T@Xr;tot=tot-cont[r]+new;cont[r]=new
        return float((tot-G).square().sum().sqrt())
z=np.load(a.init);u=torch.tensor(z['U'],**d).reshape(R,B,n);w=torch.tensor(z['w'],**d).reshape(R,B,C)
t0=time.time();print(json.dumps(dict(projected=project(u,w,a.cycles+2))),flush=True)
u.requires_grad_(True);w.requires_grad_(True);best=1e9
def loss():
    X=torch.bmm(w.transpose(1,2),feats(u)).reshape(-1,len(ii));return (X.T@X-G).square().sum()
for it in range(a.outer):
    opt=torch.optim.LBFGS([u,w],lr=1,max_iter=a.inner,history_size=50,line_search_fn='strong_wolfe')
    def closure():
        opt.zero_grad();l=loss();l.backward();return l
    opt.step(closure)
    if it%5==4:project(u.detach(),w.data,a.cycles)
    err=float(loss().detach().sqrt())
    if err<best:
        best=err;np.savez(a.out,U=(u/u.norm(dim=-1,keepdim=True)).detach().cpu().numpy(),w=w.detach().cpu().numpy(),err=err)
    print(json.dumps(dict(it=it,err=err,best=best,s=round(time.time()-t0))),flush=True)
