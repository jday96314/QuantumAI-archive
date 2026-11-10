"""Train an alternative compact-ensemble THC shape in float64.

Each ensemble is whitened independently and all physical factor columns are
optimized jointly against the full reference spectrum. Angle quantization and
normalization acceptance are separate later gates.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    p.add_argument('--steps',type=int,default=10000);p.add_argument('--lr',type=float,default=.001)
    p.add_argument('--ridge',type=float,default=1e-8);p.add_argument('--normalize-trace',action='store_true')
    p.add_argument('--optimizer',choices=['adam','lbfgs'],default='adam')
    p.add_argument('--history',type=int,default=40)
    a=p.parse_args();torch.set_num_threads(2)
    z=dict(np.load(a.input));sp=np.load('../../3/parameters/reference_spectrum.npz')
    d=dict(device='cuda',dtype=torch.float64)
    ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d)
    R,B,N=z['U'].shape;C=z['w'].shape[2]
    G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d)
    constant=G.square().sum();identity=torch.eye(B,**d)
    diagonal=torch.tensor((sp['ii']==sp['jj']).astype(float),**d)
    trace_norm=diagonal@G@diagonal
    u=torch.tensor(z['U'],requires_grad=True,**d)
    with torch.no_grad():
        pp=u[:,:,ii]*u[:,:,jj]*fac
        L=torch.linalg.cholesky(pp@pp.transpose(1,2)+1e-12*identity)
        initial=L.transpose(1,2)@torch.tensor(z['w'],**d)
    f=initial.clone().requires_grad_()

    def state():
        un=u/u.norm(dim=2,keepdim=True)
        P=un[:,:,ii]*un[:,:,jj]*fac
        L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-12*identity)
        w=torch.linalg.solve_triangular(L.transpose(1,2),f,upper=True)
        if a.normalize_trace:
            w=w*torch.sqrt(trace_norm/w.sum(1).square().sum().clamp(min=1e-20))
        X=(w.transpose(1,2)@P).reshape(R*C,-1)
        error=constant-2*(X*(X@G)).sum()+(X@X.T).square().sum()
        return error,un,w,X

    if a.optimizer=='lbfgs':
        opt=torch.optim.LBFGS([u,f],lr=a.lr,max_iter=10,max_eval=16,
                              history_size=a.history,line_search_fn='strong_wolfe',
                              tolerance_grad=1e-12,tolerance_change=1e-15)
    else:opt=torch.optim.Adam([u,f],lr=a.lr)
    best=float('inf');begin=time.time();evaluations=0
    def closure():
        nonlocal evaluations
        opt.zero_grad();error,un,w,X=state()
        objective=error+a.ridge*w.square().sum()
        objective.backward();evaluations+=1
        return objective
    for step in range(a.steps):
        if step==0:
            error,un,w,X=state()
            direct=(X.T@X-G).square().sum()
            assert abs(float((direct-error).detach()))<1e-8
        if a.optimizer=='lbfgs':opt.step(closure)
        else:
            closure();opt.step()
            opt.param_groups[0]['lr']=a.lr*(.05+.95*.5*(1+np.cos(np.pi*step/a.steps)))
        if step%(10 if a.optimizer=='lbfgs' else 100)==0 or step==a.steps-1:
            with torch.no_grad():
                error,un,w,_=state();er=float(error.clamp(min=0).sqrt())
                if er<best:
                    best=er
                    path=Path(a.output);tmp=path.with_name(path.stem+'.tmp.npz')
                    np.savez(tmp,U=un.cpu().numpy(),w=w.cpu().numpy(),W=np.zeros((R,C)),tree_edges=z['tree_edges'],err2=er)
                    tmp.replace(path)
                print(json.dumps(dict(step=step,evaluations=evaluations,optimizer=a.optimizer,error_2=er,best=best,R=R,B=B,C=C,weight_squared_norm=float(w.square().sum()),seconds=time.time()-begin)),flush=True)
    print(json.dumps(dict(status='FINISHED',scope='unquantized alternative ensemble model; no native acceptance',best_error_2=best,seconds=time.time()-begin)),flush=True)


if __name__=='__main__':main()
