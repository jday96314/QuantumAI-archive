"""Bounded float64 physical tensor fit of multiblock overcomplete dictionaries.

Whiten each *symmetric feature* Gram before parameterizing its covariance.
The full target spectrum, including negative modes, remains in the loss.
"""
import argparse
import time
import json
from pathlib import Path
import numpy as np
import torch
from common import HERE,SPECTRUM,report

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--steps',type=int,default=2000);p.add_argument('--lr',type=float,default=.002)
    p.add_argument('--optimizer',choices=['adam','lbfgs'],default='adam');p.add_argument('--ridge',type=float,default=1e-10)
    p.add_argument('--trace',action='store_true');p.add_argument('--history',type=int,default=30)
    a=p.parse_args();torch.set_num_threads(2);begin=time.time()
    z=dict(np.load(a.input));sp=np.load(SPECTRUM);R,B,N=z['U'].shape;C=z['w'].shape[2]
    d=dict(device='cuda',dtype=torch.float64);ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d)
    G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d);constant=G.square().sum()
    eye=torch.eye(B,**d);diag=torch.tensor((ii==jj).astype(float),**d);trace_target=diag@G@diag
    u=torch.tensor(z['U'],requires_grad=True,**d)
    with torch.no_grad():
        P=u[:,:,ii]*u[:,:,jj]*fac
        L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*eye)
        f0=L.transpose(1,2)@torch.tensor(z['w'],**d)
    f=f0.clone().requires_grad_()
    def state():
        un=u/u.norm(dim=2,keepdim=True)
        P=un[:,:,ii]*un[:,:,jj]*fac
        L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*eye)
        w=torch.linalg.solve_triangular(L.transpose(1,2),f,upper=True)
        if a.trace:w=w*torch.sqrt(trace_target/w.sum(1).square().sum().clamp(min=1e-20))
        X=(w.transpose(1,2)@P).reshape(R*C,-1)
        error=constant-2*(X*(X@G)).sum()+(X@X.T).square().sum()
        return error,un,w,X
    if a.optimizer=='adam':opt=torch.optim.Adam([u,f],lr=a.lr)
    else:opt=torch.optim.LBFGS([u,f],lr=a.lr,max_iter=5,max_eval=8,history_size=a.history,line_search_fn='strong_wolfe',tolerance_grad=1e-12,tolerance_change=1e-15)
    evals=0;best=float('inf');records=[]
    def closure():
        nonlocal evals
        opt.zero_grad();err,un,w,X=state();loss=err+a.ridge*w.square().sum();loss.backward();evals+=1
        return loss
    with torch.no_grad():
        er,_,_,X=state();direct=(X.T@X-G).square().sum()
        initial=float(direct.sqrt());initial_discrepancy=abs(float(er-direct))
        assert initial_discrepancy<2e-9
    for step in range(a.steps):
        if a.optimizer=='adam':
            closure();opt.step();opt.param_groups[0]['lr']=a.lr*(.05+.95*.5*(1+np.cos(np.pi*step/a.steps)))
        else:opt.step(closure)
        if step%50==0 or step==a.steps-1:
            with torch.no_grad():
                er,un,w,X=state();direct=(X.T@X-G).norm();err=float(direct)
                discrepancy=abs(float(er-direct.square()))
                assert discrepancy<5e-9,discrepancy
                if err<best:
                    best=err
                    out=dict(z,U=un.cpu().numpy(),w=w.cpu().numpy(),err2=err)
                    out.pop('angle_words',None)
                    path=Path(a.output);tmp=path.with_name(path.stem+'.tmp.npz');np.savez(tmp,**out);tmp.replace(path)
                rec=dict(step=step,evaluations=evals,error_2=err,best=best,weight_squared_norm=float(w.square().sum()),seconds=time.time()-begin)
                records.append(rec);print(json.dumps(rec),flush=True)
    report(HERE/f'results/{Path(a.output).stem}_relax.json',dict(status='CONTINUOUS_FIT_ONLY',input=a.input,output=a.output,
        R=R,B=B,C=C,optimizer=a.optimizer,steps=a.steps,ridge=a.ridge,trace_normalized=a.trace,
        initial_error_2=initial,initial_direct_contraction_squared_discrepancy=initial_discrepancy,best_error_2=best,records=records,seconds=time.time()-begin))

if __name__=='__main__':main()
