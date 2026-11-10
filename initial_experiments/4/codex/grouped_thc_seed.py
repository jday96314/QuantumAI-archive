"""Initialize small THC ensembles from physically similar reference groups.

Agglomerate by covariance overlap, retain each group's principal physical
factor columns, and select its shared rank-one features by simultaneous OMP.
This is an alternative factorization family, not a submitted circuit.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--groups',type=int,default=10)
    p.add_argument('--basis',type=int,default=64)
    p.add_argument('--columns',type=int,default=35)
    p.add_argument('--output',default='parameters/grouped10_b64_c35_seed.npz')
    a=p.parse_args();torch.set_num_threads(2);begin=time.time()
    ref=np.load('../../3/parameters/chem_reference.npz')
    sp=np.load('../../3/parameters/reference_spectrum.npz')
    d=dict(device='cuda',dtype=torch.float64)
    u=torch.tensor(ref['U'],**d);w=torch.tensor(ref['w'],**d)
    ii,jj=sp['ii'],sp['jj'];fac=torch.tensor(sp['fac'],**d)
    P=u[:,:,ii]*u[:,:,jj]*fac
    X=w.transpose(1,2)@P
    oldR,C,D=X.shape
    cross=(X.reshape(-1,D)@X.reshape(-1,D).T).reshape(oldR,C,oldR,C)
    overlap=cross.square().sum((1,3)).cpu().numpy()
    groups=[[i] for i in range(oldR)];maximum=int(np.ceil(oldR/a.groups))
    while len(groups)>a.groups:
        best=None
        for i in range(len(groups)):
            for j in range(i):
                left,right=groups[i],groups[j]
                if len(left)+len(right)>maximum:continue
                denom=np.sqrt(overlap[np.ix_(left,left)].sum()*overlap[np.ix_(right,right)].sum())
                score=overlap[np.ix_(left,right)].sum()/max(denom,1e-30)
                if best is None or score>best[0]:best=(score,i,j)
        if best is None:
            maximum+=1
            continue
        _,i,j=best;groups[j]+=groups[i];groups.pop(i)
    outs=[];weights=[];records=[]
    for members in groups:
        target=X[members].reshape(-1,D)
        _,v=torch.linalg.eigh(target@target.T)
        target=v[:,-a.columns:].T@target
        candidates=u[members].reshape(-1,76)
        if len(candidates)<a.basis:
            # A singleton reference ensemble has only60 atoms. Extra
            # independent candidates permit a64-atom output without cloning
            # a feature and making its whitening Gram singular.
            extra=np.load('candidates/norm_cycle2_signed/submission.npz')['U'].reshape(-1,76)
            candidates=torch.cat([candidates,torch.tensor(extra,**d)])
        pc=candidates[:,ii]*candidates[:,jj]*fac
        gram=pc@pc.T;initial_gram=gram.clone();rhs=pc@target.T
        selected=[]
        for _ in range(a.basis):
            diag=gram.diag().clamp(min=1e-16)
            gain=rhs.square().sum(1)/diag
            gain[diag<1e-10]=-1
            if selected:gain[selected]=-1
            j=int(gain.argmax());assert float(gain[j])>=0
            selected.append(j)
            vj=gram[:,j]/diag[j].sqrt();hj=rhs[j]/diag[j].sqrt()
            gram=gram-torch.outer(vj,vj);rhs=rhs-torch.outer(vj,hj)
        inds=torch.tensor(selected,device='cuda')
        features=pc[inds]
        ww=torch.linalg.solve(initial_gram[inds][:,inds]+1e-12*torch.eye(a.basis,**d),features@target.T)
        residual=ww.T@features-target
        records.append(dict(members=members,selected_atoms=selected,target_projection_error=float(residual.square().sum().sqrt())))
        outs.append(candidates[inds].cpu().numpy());weights.append(ww.cpu().numpy())
    U=np.array(outs);W=np.array(weights)
    model=torch.tensor(W,**d).transpose(1,2)@(torch.tensor(U[:,:,ii]*U[:,:,jj]*sp['fac'],**d))
    model=model.reshape(a.groups*a.columns,D)
    G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d)
    error=(G.square().sum()-2*(model*(model@G)).sum()+(model@model.T).square().sum()).clamp(min=0).sqrt()
    tree=np.load('candidates/norm_cycle2_signed/submission.npz')['tree_edges']
    np.savez(a.output,U=U,w=W,W=np.zeros((a.groups,a.columns)),tree_edges=tree,err2=float(error))
    report=dict(status='INITIALIZED',scope='classical unquantized alternative ensemble factorization',R=a.groups,B=a.basis,C=a.columns,error_2=float(error),groups=records,seconds=time.time()-begin)
    Path('results/'+Path(a.output).stem+'.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='groups'}),flush=True)


if __name__=='__main__':main()
