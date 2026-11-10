"""Fixed-word gradient refitting along physical constraint tangents.

Use the actual residual gradient, then project the direction against the
linearized gap and normalization limits, retaining their feasible headroom.
Every accepted step is independently assessed with float64 decoded physics.
"""
import argparse,json,time
from pathlib import Path
import torch
import numpy as np
from audit_core import Metric,load,save
from chem_rank_screen import HERE,quantities


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True);p.add_argument('--steps',type=int,default=500)
    p.add_argument('--lambda-cap',type=float,default=320.);p.add_argument('--lambda-eff-cap',type=float,default=51.);p.add_argument('--inward',type=float,default=.002);p.add_argument('--tag',default='tangent_refit')
    p.add_argument('--ridge',type=float,default=0.);p.add_argument('--step-norm',type=float,default=.05);p.add_argument('--project-constraints',action='store_true');args=p.parse_args()
    torch.set_num_threads(1);z,u,w,tree,words=load(args.init);m=Metric(tree)
    t=m.tensor(words.reshape(1,-1,75));u=m.decode(t);weights=m.tensor(w);offsets=m.tensor(z['W'])
    m.project_gap(u,weights,offsets,target=3.91515);history=[];started=time.time()
    chol=None
    if args.ridge:
        feature=m.feats(u)[0];gram_feature=feature@feature.T
        chol=torch.linalg.cholesky(gram_feature+args.ridge*torch.eye(len(feature),**m.d))
    def precondition(values):
        return [torch.cholesky_solve(values[0][0],chol)[None],values[1]] if chol is not None else list(values)
    def dot(a,b):return sum((x*y).sum() for x,y in zip(a,b))
    def project_constraints(ww,oo):
        for _ in range(8):
            with torch.enable_grad():
                ww=ww.detach().requires_grad_(True);oo=oo.detach().requires_grad_(True)
                g,l,e=m.physical(u,ww,oo)
                violations=[(3.91510-g)/3.91510,(l-args.lambda_cap+1e-6)/args.lambda_cap,
                            (e-args.lambda_eff_cap**2+1e-5)/args.lambda_eff_cap**2]
                index=max(range(3),key=lambda i:float(violations[i].detach()))
                value=float(violations[index].detach())
                if value<=0:return ww.detach(),oo.detach()
                normal=torch.autograd.grad(violations[index],(ww,oo))
            scale=value/max(float(dot(normal,normal)),1e-30)
            with torch.no_grad():
                old=(ww.detach(),oo.detach());changed=False
                for damping in (1.,.5,.25,.125):
                    nw,no=[v-damping*scale*n for v,n in zip(old,normal)]
                    ng,nl,ne=m.physical(u,nw,no)
                    merit=max(float((3.91510-ng)/3.91510),float((nl-args.lambda_cap+1e-6)/args.lambda_cap),
                              float((ne-args.lambda_eff_cap**2+1e-5)/args.lambda_eff_cap**2))
                    if merit<value:ww,oo=nw,no;changed=True;break
                if not changed:return old
        return ww.detach(),oo.detach()
    for step in range(args.steps+1):
        current=quantities(m,u,weights,offsets);row=dict(step=step,seconds=time.time()-started,**current)
        history.append(row);print(json.dumps(row),flush=True)
        save(HERE/'parameters'/f'{args.tag}_last.npz',u,weights,tree,t,W=offsets.cpu().numpy())
        (HERE/'results'/f'{args.tag}_history.json').write_text(json.dumps(history,indent=2)+'\n')
        if step==args.steps:break
        ww=weights.detach().requires_grad_(True);oo=offsets.detach().requires_grad_(True)
        gradient=torch.autograd.grad(m.loss(u,ww),ww)[0]
        gap,lam,effective=m.physical(u,ww,oo)
        normal_gap=torch.autograd.grad(gap,(ww,oo),retain_graph=True)
        normal_lam=torch.autograd.grad(lam,(ww,oo),retain_graph=True)
        normal_eff=torch.autograd.grad(effective,(ww,oo))
        with torch.no_grad():
            direction=precondition([gradient,torch.zeros_like(offsets)])
            factor=-args.step_norm/dot(direction,direction).sqrt().clamp(min=1e-30)
            direction=[factor*v for v in direction]
            # Project onto all three linearized halfspaces. Enumerating the
            # eight possible active sets solves this small convex QP exactly.
            normals=[[-v for v in normal_gap],list(normal_lam),list(normal_eff)]
            bounds=[current['gap']-3.91515,args.lambda_cap-current['Lambda'],args.lambda_eff_cap**2-current['lambda_eff']**2]
            if bounds[1]<.01:bounds[1]=-args.inward
            if bounds[2]<.1:bounds[2]=-5*args.inward
            prec_normals=[precondition(v) for v in normals]
            norms=[float(dot(v,pv).clamp(min=1e-30).sqrt()) for v,pv in zip(normals,prec_normals)]
            normals=[[x/max(norm,1e-30) for x in values] for norm,values in zip(norms,normals)]
            prec_normals=[[x/max(norm,1e-30) for x in values] for norm,values in zip(norms,prec_normals)]
            bounds=np.asarray([v/max(norm,1e-30) for v,norm in zip(bounds,norms)])
            gram=np.asarray([[float(dot(a,b)) for b in prec_normals] for a in normals])
            old=np.asarray([float(dot(v,direction)) for v in normals]);best=None
            for active_bits in range(8):
                active=[i for i in range(3) if active_bits>>i&1];multipliers=np.zeros(3)
                if active:
                    multipliers[active]=np.linalg.lstsq(gram[np.ix_(active,active)],(old-bounds)[active],rcond=1e-12)[0]
                if multipliers.min()<-1e-10:continue
                if np.any(old-gram@multipliers>bounds+1e-9):continue
                cost=float(multipliers@gram@multipliers)
                if best is None or cost<best[0]:best=(cost,multipliers)
            if best is None:
                print(json.dumps(dict(status='NO_LINEARIZED_FEASIBLE_DIRECTION',step=step)),flush=True);break
            direction=[v-sum(float(multiplier)*normal[j] for multiplier,normal in zip(best[1],prec_normals)) for j,v in enumerate(direction)]
            accepted=False
            for exponent in range(18):
                rate=2.**(-exponent);trial_w=weights+rate*direction[0];trial_o=offsets+rate*direction[1]
                error=float(m.loss(u,trial_w).sqrt())
                if error>=current['error2']-1e-11:continue
                m.project_gap(u,trial_w,trial_o,target=3.91515)
                if args.project_constraints:
                    trial_w,trial_o=project_constraints(trial_w,trial_o)
                    if float(m.loss(u,trial_w).sqrt())>=current['error2']-1e-11:continue
                tg,tl,te=m.physical(u,trial_w,trial_o)
                if float(tg)<3.91504 or float(tl)>args.lambda_cap or float(te.sqrt())>args.lambda_eff_cap:continue
                weights,offsets=trial_w,trial_o;accepted=True
                save(HERE/'parameters'/f'{args.tag}_best.npz',u,weights,tree,t,W=offsets.cpu().numpy())
                break
            if not accepted:
                print(json.dumps(dict(status='NO_FEASIBLE_DESCENT',step=step,direction_norm=float(direction[0].norm()),gap_direction=float(dot(normal_gap,direction)),lambda_direction=float(dot(normal_lam,direction)))),flush=True)
                break
