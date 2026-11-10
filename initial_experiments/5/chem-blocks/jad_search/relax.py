"""Bounded physical Gram relaxation of an overcomplete multiblock dictionary.

Uses per-ensemble projector-Gram whitening, never Euclidean PSD projection in
unwhitened atom coordinates. Full signed reference spectrum is retained.
"""
import argparse,time
import numpy as np
import torch
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--steps',type=int,default=300);p.add_argument('--lr',type=float,default=.5)
    p.add_argument('--max-iter',type=int,default=5);p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--mode',choices=['tensor','factor'],default='tensor');p.add_argument('--ridge',type=float,default=1e-10)
    a=p.parse_args();torch.set_num_threads(2);beg=time.time();z=dict(np.load(a.input));_,sp=load()
    R,B,N=z['U'].shape;C=z['w'].shape[2];d=dict(device='cuda',dtype=torch.float64)
    ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d);G=torch.tensor(gram(sp),**d);constant=G.square().sum();eye=torch.eye(B,**d)
    u=torch.tensor(z['U'],requires_grad=True,**d)
    with torch.no_grad():
        P=u[:,:,ii]*u[:,:,jj]*fac;L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*eye)
        init=L.transpose(1,2)@torch.tensor(z['w'],**d)
    f=init.clone().requires_grad_();target=torch.tensor(z.get('target',np.zeros((R,C,len(ii)))),**d)
    def state():
        U=u/u.norm(dim=2,keepdim=True);P=U[:,:,ii]*U[:,:,jj]*fac
        L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*eye)
        w=torch.linalg.solve_triangular(L.transpose(1,2),f,upper=True)
        X=(w.transpose(1,2)@P).reshape(R*C,-1)
        err=constant-2*(X*(X@G)).sum()+(X@X.T).square().sum()
        local=(X.reshape(R,C,-1)-target).square().sum()
        return err,local,U,w,X
    opt=torch.optim.LBFGS([u,f],lr=a.lr,max_iter=a.max_iter,max_eval=2*a.max_iter,
        history_size=30,line_search_fn='strong_wolfe',tolerance_grad=1e-12,tolerance_change=1e-15)
    neval=0;best=float('inf');records=[]
    def closure():
        nonlocal neval
        opt.zero_grad();err,local,U,w,X=state()
        loss=(err if a.mode=='tensor' else local)+a.ridge*w.square().sum()
        loss.backward();neval+=1;return loss
    with torch.no_grad():
        err,_,_,_,X=state();direct=(X.T@X-G).square().sum()
        assert abs(float(direct-err))<1e-8
        initial=float(direct.sqrt())
    for step in range(a.steps):
        opt.step(closure)
        if step%10==0 or step==a.steps-1 or time.time()-beg>a.seconds:
            with torch.no_grad():
                err,local,U,w,X=state();er=float(err.clamp(min=0).sqrt())
                if er<best:
                    best=er
                    np.savez(a.output,U=U.cpu().numpy(),w=w.cpu().numpy(),offset=z['offset'],target=z.get('target',np.zeros((R,C,len(ii)))),tree_edges=z['tree_edges'],err2=er)
                row=dict(step=step,evaluations=neval,error_2=er,best_error_2=best,factor_projection_error=float(local.sqrt()),weight_squared_norm=float(w.square().sum()),seconds=time.time()-beg)
                records.append(row);emit(row)
            if time.time()-beg>a.seconds:break
    zz=np.load(a.output);Xn=np.einsum('rbc,rbd->rcd',zz['w'],feature(zz['U'],sp)).reshape(R*C,-1)
    direct=float(np.linalg.norm(Xn.T@Xn-gram(sp)))
    assert abs(direct-best)<1e-6
    data=dict(status='CONTINUOUS_RESEARCH_ONLY',shape=[R,B,C],input=a.input,output=a.output,initial_error_2=initial,best_error_2=direct,mode=a.mode,records=records,seconds=time.time()-beg)
    report(Path(a.output).stem+'_relax.json',data);emit(data)

if __name__=='__main__':main()
