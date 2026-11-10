"""Exact-word pair descent ranked by the coupled Gauss-Newton matrix.

All 75 choose 2 angle pairs and four signs are screened per row. Exact
single-row tensor changes are then evaluated by 76x76 quadratic forms,
followed by full float64 native-metric and physical checks for accepted
batches. Unlike independent gradient proposals, both signs of a pair are
eligible even when either individual move is uphill. No verifier changes.
"""
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric,load,save
from chem_rank_screen import quantities
HERE=Path(__file__).resolve().parent

def jacobian(m,t,tree):
 B=len(t);v=torch.zeros((B,76),**m.d);v[:,0]=1;J=torch.zeros((B,75,76),**m.d)
 for j,(p,c) in enumerate(tree):
  co=(t[:,j]*np.pi/128).cos();sn=(t[:,j]*np.pi/128).sin()
  jp,jc=J[:,:,p].clone(),J[:,:,c].clone();J[:,:,p]=co[:,None]*jp-sn[:,None]*jc;J[:,:,c]=sn[:,None]*jp+co[:,None]*jc
  vp,vc=v[:,p].clone(),v[:,c].clone();v[:,p]=co*vp-sn*vc;v[:,c]=sn*vp+co*vc
  J[:,j,p]=-np.pi/128*v[:,c];J[:,j,c]=np.pi/128*v[:,p]
 return v,J

@torch.no_grad()
def main():
 p=argparse.ArgumentParser();p.add_argument('--init',type=Path,required=True);p.add_argument('--steps',type=int,default=100);p.add_argument('--tag',required=True)
 p.add_argument('--pairs',type=int,default=12);p.add_argument('--lambda-cap',type=float,default=330);p.add_argument('--lambda-eff-cap',type=float,default=51)
 p.add_argument('--refit',action='store_true');a=p.parse_args();torch.set_num_threads(1);started=time.time()
 z,u0,w,tree,words=load(a.init);assert w.shape[0]==1;m=Metric(tree);t=m.tensor(words.reshape(1,-1,75));w=m.tensor(w);off=m.tensor(z['W']);B=w.shape[1]
 pair=torch.triu_indices(75,75,1,device=t.device);pi,pj=pair;sgn=t.new_tensor([[1,1],[1,-1],[-1,1],[-1,-1]])
 rows=[];accepted=0;best=float('inf');refits=0;pred_discrepancies=[]
 def feasible(uu,ww,oo):
  g,l,e=m.physical(uu,ww,oo)
  return float(l)<=a.lambda_cap+1e-8 and float(g)>=3.91504 and float(e.sqrt())<=a.lambda_eff_cap+1e-8
 for step in range(a.steps+1):
  u=m.decode(t)
  if a.refit and step:
   wp,_=m.project(u,w,1)
   # Preserve physical normalization while taking a PSD weight refit.
   for scale in (1.,.5,.25,.125,.0625,.015625,.00390625):
    wt=w+scale*(wp-w);ot=off.clone();m.project_gap(u,wt,ot)
    if feasible(u,wt,ot):w,off=wt,ot;refits+=1;break
  metrics=quantities(m,u,w,off);err=metrics['error2'];rec=dict(step=step,accepted_rows=accepted,weight_refits=refits,**metrics,seconds=time.time()-started)
  rows.append(rec);print(json.dumps(rec),flush=True)
  save(HERE/'parameters'/f'{a.tag}_last.npz',u,w,tree,t,W=off.cpu().numpy())
  if err<best:
   best=err;save(HERE/'parameters'/f'{a.tag}_best_error.npz',u,w,tree,t,W=off.cpu().numpy())
  out=dict(status='FINITE_COUPLED_WORD_SCREEN',arguments={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,max_prediction_discrepancy=max(pred_discrepancies,default=0),finished=step==a.steps)
  (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(out,indent=2)+'\n')
  if step==a.steps:break
  P=m.feats(u)[0];X=w[0].T@P;D=X.T@X-m.G;Z=w[0]@X
  L=m.unpack(Z@D);Zm=m.unpack(Z);alpha=w[0].square().sum(1)
  uv,J=jacobian(m,t[0],tree);Lu=(L@uv[:,:,None]).squeeze(-1);Zu=(Zm@uv[:,:,None]).squeeze(-1)
  grad=8*(J@Lu[:,:,None]).squeeze(-1);v=(J@Zu[:,:,None]).squeeze(-1)
  Hess=8*Z.square().sum(1)[:,None,None]*(J@J.transpose(-1,-2))+16*v[:,:,None]*v[:,None,:]
  hd=Hess.diagonal(dim1=-2,dim2=-1);base=.5*(hd[:,pi]+hd[:,pj]);cross=Hess[:,pi,pj]
  score=base[:,:,None]+grad[:,pi,None]*sgn[:,0]+grad[:,pj,None]*sgn[:,1]+cross[:,:,None]*(sgn[:,0]*sgn[:,1])
  vals,ix=score.flatten(1).topk(a.pairs,dim=1,largest=False);ip=ix//4;isign=ix%4
  delta=torch.zeros((B,a.pairs+8,75),**m.d);batch=torch.arange(B,device=t.device)[:,None];cn=torch.arange(a.pairs,device=t.device)[None,:]
  delta[batch,cn,pi[ip]]=sgn[isign,0];delta[batch,cn,pj[ip]]=sgn[isign,1]
  sv,si=(.5*hd[:,:,None]+grad[:,:,None]*t.new_tensor([1,-1])).flatten(1).topk(8,dim=1,largest=False)
  delta[batch,torch.arange(8,device=t.device)[None,:]+a.pairs,si//2]=1-2*(si%2).double()
  candidates=t[0,:,None,:]+delta;cu=m.decode(candidates);CP=m.feats(cu);dp=CP-P[:,None,:]
  # Exact row-local residual update, including every quadratic interaction.
  norm=dp.square().sum(-1);dot=(dp*Z[:,None,:]).sum(-1)
  change=4*(dp*(Z@D)[:,None,:]).sum(-1)+2*Z.square().sum(1)[:,None]*norm+2*dot.square()+4*alpha[:,None]*dot*norm+alpha[:,None].square()*norm.square()
  values,choice=change.min(1);rank=values.argsort();moves=delta[torch.arange(B,device=t.device),choice]
  rec['minimum_pair_GN']=float(vals.min());rec['minimum_exact_row_change']=float(values.min());rec['predicted_improving_rows']=int((values<-1e-13).sum())
  changed=False
  groups=[rank[:min(k,B)] for k in (128,64,32,16,8,4,2,1)];groups += [rank[i:i+1] for i in range(1,16)]
  for group in groups:
   trial=t.clone();trial[0,group]+=moves[group];uu=m.decode(trial);ne2=float(m.loss(uu,w));expected=err**2+float(values[group].sum())
   if len(group)==1:
    discrepancy=abs(ne2-expected);pred_discrepancies.append(discrepancy);assert discrepancy<1e-9,discrepancy
   if ne2>=err**2-1e-13:continue
   oo=off.clone();m.project_gap(uu,w,oo)
   if not feasible(uu,w,oo):continue
   t,off=trial,oo;accepted+=len(group);changed=True;break
  if not changed:
   rec['status']='NO_ACCEPTED_PAIR_BATCH';out['finished']=True;out['termination']='NO_ACCEPTED_PAIR_BATCH';out['max_prediction_discrepancy']=max(pred_discrepancies,default=0)
   (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rec),flush=True);break
 print(json.dumps(dict(status='COMPLETE_FINITE_PAIR_SEARCH',best_error=best,seconds=time.time()-started)),flush=True)
if __name__=='__main__':main()
