"""Exact integer row descent while preserving a low-normalization weight gauge.

Float32 derivatives only rank proposals. Every accepted dictionary is decoded
and assessed in float64 with the actual signed tensor and physical formulas.
Unlike the prior refiner, word trials need not accept an unconstrained PSD
weight refit, which can destroy their otherwise feasible normalization.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
import torch
from audit_core import Metric,load,save
from chem_rank_screen import HERE,quantities


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--init',type=Path,required=True)
    parser.add_argument('--steps',type=int,default=100)
    parser.add_argument('--tag',default='fixed_weight_repair')
    parser.add_argument('--lambda-cap',type=float,default=320.)
    parser.add_argument('--lambda-eff-cap',type=float,default=51.)
    parser.add_argument('--refit',action='store_true')
    parser.add_argument('--refit-margin',type=float,default=0.)
    parser.add_argument('--normalization-penalty',type=float,default=0.)
    args=parser.parse_args();torch.set_num_threads(1);torch.manual_seed(20261025)
    z,u,w,tree,words=load(args.init);assert w.shape[0]==1
    exact=Metric(tree);fast=copy.copy(exact);fast.d=dict(device='cuda',dtype=torch.float32)
    for name in ('fac','G','h0'):setattr(fast,name,getattr(exact,name).float())
    t=exact.tensor(words.reshape(1,-1,75));weights=exact.tensor(w);offsets=exact.tensor(z['W'])
    started=time.time();history=[];accepted=0;best_error=float('inf')
    def feasible(u,w,o,margin=0.):
        gap,lam,effective=exact.physical(u,w,o)
        return float(lam)<=args.lambda_cap-margin and float(gap)>=3.91504 and float(effective.sqrt())<=args.lambda_eff_cap
    for step in range(args.steps+1):
        u=exact.decode(t)
        if args.refit:
            with torch.no_grad():
                proposal,_=exact.project(u,weights,1)
                lo,hi=0.,1.
                for _ in range(24):
                    a=(lo+hi)/2;trial=weights+a*(proposal-weights)
                    if feasible(u,trial,offsets,args.refit_margin):lo=a
                    else:hi=a
                weights=weights+lo*(proposal-weights)
        with torch.no_grad():
            err=float(exact.loss(u,weights).sqrt());metrics=quantities(exact,u,weights,offsets)
            record=dict(step=step,accepted_rows=accepted,seconds=time.time()-started,**metrics)
            objective=err**2+args.normalization_penalty*(metrics['lambda_eff']/args.lambda_eff_cap)**2
            record['objective']=objective
            history.append(record);print(json.dumps(record),flush=True)
            save(HERE/'parameters'/f'{args.tag}_last.npz',u,weights,tree,t,W=offsets.cpu().numpy())
            if err<best_error:
                best_error=err
                save(HERE/'parameters'/f'{args.tag}_best_error.npz',u,weights,tree,t,W=offsets.cpu().numpy())
            (HERE/'results'/f'{args.tag}_history.json').write_text(json.dumps(history,indent=2)+'\n')
        if step==args.steps:break
        tf=t.float().detach().requires_grad_(True);wf=weights.float()
        uf=fast.decode(tf);loss=fast.loss(uf,wf)
        gradient=torch.autograd.grad(loss,tf,retain_graph=bool(args.normalization_penalty))[0].detach()[0]
        norm_gradient=None
        if args.normalization_penalty:
            effective=fast.physical(uf,wf,offsets.float())[2]/args.lambda_eff_cap**2
            norm_gradient=torch.autograd.grad(effective,tf)[0].detach()[0]
            gradient=gradient+args.normalization_penalty*norm_gradient
        with torch.no_grad():
            P=fast.feats(uf)[0];X=wf[0].T@P;D=X.T@X-fast.G;Z=(wf[0]@wf[0].T)@P
            tf=tf.detach()[0];B=len(tf)
            v=torch.zeros((B,76),**fast.d);v[:,0]=1
            J=torch.zeros((B,75,76),**fast.d)
            for j,(p,c) in enumerate(tree):
                co=(tf[:,j]*np.pi/128).cos();sn=(tf[:,j]*np.pi/128).sin()
                jp,jc=J[:,:,p].clone(),J[:,:,c].clone()
                J[:,:,p]=co[:,None]*jp-sn[:,None]*jc;J[:,:,c]=sn[:,None]*jp+co[:,None]*jc
                vp,vc=v[:,p].clone(),v[:,c].clone();v[:,p]=co*vp-sn*vc;v[:,c]=sn*vp+co*vc
                J[:,j,p]=-np.pi/128*v[:,c];J[:,j,c]=np.pi/128*v[:,p]
            h=(8*J.square().sum(-1)*Z.square().sum(-1)[:,None]).clamp(min=1e-18)
            order=(.5*h-gradient.abs()).argsort(-1);direction=-gradient.sign()
            candidates=[]
            for count in (1,2,4,8):
                ix=order[:,:count];delta=torch.zeros_like(tf)
                delta.scatter_(1,ix,direction.gather(1,ix));candidates.append(tf+delta)
            candidates=torch.stack(candidates,1)
            dp=fast.feats(fast.decode(candidates))-P[:,None,:]
            norm=dp.square().sum(-1);dot=(dp*Z[:,None,:]).sum(-1)
            diagonal=wf[0].square().sum(-1)[:,None]
            prediction=4*(dp*(Z@D)[:,None,:]).sum(-1)+2*Z.square().sum(-1)[:,None]*norm+2*dot.square()+4*diagonal*dot*norm+diagonal.square()*norm.square()
            if norm_gradient is not None:prediction+=args.normalization_penalty*((candidates-tf[:,None,:])*norm_gradient[:,None,:]).sum(-1)
            prediction[norm<1e-20]=torch.inf
            values,choice=prediction.min(-1);ranking=values.argsort()
            moves=candidates[torch.arange(B,device=tf.device),choice]-tf
            changed=False
            groups=[ranking[:min(count,B)] for count in (128,64,32,16,8,4,2,1)]
            groups += [ranking[j:j+1] for j in range(1,16)]
            for group in groups:
                trial=t.clone();trial[0,group]+=moves[group].double()
                uu=exact.decode(trial);ne=float(exact.loss(uu,weights).sqrt())
                if not args.normalization_penalty and ne>=err-1e-10:continue
                off=offsets.clone();exact.project_gap(uu,weights,off,target=3.91515)
                if not feasible(uu,weights,off):continue
                if args.normalization_penalty:
                    effective=float(exact.physical(uu,weights,off)[2])
                    if ne**2+args.normalization_penalty*effective/args.lambda_eff_cap**2>=objective-1e-12:continue
                t=trial;offsets=off;accepted+=len(group);changed=True
                save(HERE/'parameters'/f'{args.tag}_best.npz',uu,weights,tree,t,W=offsets.cpu().numpy())
                break
            if not changed:
                print(json.dumps(dict(status='NO_ACCEPTED_PROPOSAL',step=step,best_predicted=float(values.min()),candidates_tested=len(groups))),flush=True)
                break
