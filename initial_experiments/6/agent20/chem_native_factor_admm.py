"""Finite ADMM between a fixed native tensor factor and exact angle words.

The factor projection is a global orthogonal Procrustes rotation, preserving
the selected tensor. The dictionary projection uses exact integer words,
closed-form factor least squares and nonlinear row moves. This differs from
direct tensor-norm descent: ADMM dual variables can leave its local basin.
No circuit or physical acceptance is asserted by this research screen.
"""
import argparse,json,time,hashlib
from pathlib import Path
from audit_core import *
HERE=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--init',default='parameters/positive_native_integer_R4B240.npz');p.add_argument('--factor',default='parameters/psd_generator_admm_native_C352.npz');p.add_argument('--steps',type=int,default=100);p.add_argument('--seed',type=int,default=2026102620);p.add_argument('--metric-power',type=float,default=0.);p.add_argument('--metric-ridge',type=float,default=.005);p.add_argument('--dual-step',type=float,default=1.);p.add_argument('--tag',default='native_factor_admm_R4B240');p.add_argument('--global-atoms',action='store_true');p.add_argument('--inner-steps',type=int,default=1);a=p.parse_args();torch.set_num_threads(1);torch.manual_seed(a.seed);started=time.time()
 z,u,w,tree,words=load(a.init);m=Metric(tree);R,B,C=w.shape;fz=dict(np.load(a.factor));X=m.tensor(fz['factor']);assert X.shape[0]==R*C
 if words is None:words=encode(u.reshape(-1,76),tree,8,True)
 t=m.tensor(np.where(words>=128,words-256,words).reshape(R,B,75));wt=m.tensor(w);dual=torch.zeros_like(X);off=m.tensor(z['W']);rows=[];best=float('inf');accepted=0
 if a.metric_power:
  with torch.no_grad():
   _,sv,V=torch.linalg.svd(X,full_matrices=False);V=V.T;metric_weights=sv.pow(2*a.metric_power)
  def metric(T):return a.metric_ridge*T+((T@V)*metric_weights)@V.T
  mean_metric=a.metric_ridge+float(metric_weights.sum())/X.shape[1]
  metric_lipschitz=a.metric_ridge+float(metric_weights.max())
 else:
  def metric(T):return T
  mean_metric=1.
  metric_lipschitz=1.
 metric_X=metric(X)
 def fit(tt,target):
  uu=m.decode(tt);P=m.feats(uu);MP=metric(P);S=MP@P.transpose(-1,-2);ridge=1e-12*torch.eye(B,**m.d);L=torch.linalg.cholesky(S+ridge);ww=torch.cholesky_solve(MP@target.reshape(R,C,-1).transpose(1,2),L);Z=ww.transpose(1,2)@P;res=Z-target.reshape(R,C,-1)
  return uu,P,ww,Z.flatten(0,1),float((res*metric(res)).sum())
 def record(step,Y,Z,loss,changed):
  nonlocal best
  with torch.no_grad():
   err=float((Z.T@Z-m.G).norm());u=m.decode(t)
   cp,_=checkpoint(m,u,wt,off,tree,1)
   row=dict(step=step,tensor_error=err,factor_residual=float((Y-Z).norm()),dual_norm=float(dual.norm()),projection_loss=loss,changed=changed,accepted_rows=accepted,weight_norm=float(wt.norm()),metrics=cp,seconds=time.time()-started);rows.append(row);print(json.dumps(row),flush=True)
   if err<best:
    best=err;save(HERE/'parameters'/f'{a.tag}_best.npz',u,wt,tree,t,err2=err,W=off.cpu().numpy())
   save(HERE/'parameters'/f'{a.tag}_last.npz',u,wt,tree,t,err2=err,W=off.cpu().numpy())
   (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(dict(status='FINITE_NATIVE_FACTOR_ADMM',arguments=vars(a),rows=rows,best_tensor_error=best,finished=step==a.steps,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='exact-word dictionary, signed factor coefficients, physical metrics reported; no native circuit'),indent=2)+'\n')
 with torch.no_grad():
  Z=(wt.transpose(1,2)@m.feats(m.decode(t))).flatten(0,1);left,sv,right=torch.linalg.svd(Z@metric_X.T,full_matrices=False);Y=(left@right)@X;record(-1,Y,Z,float(((Y-Z)*metric(Y-Z)).sum()),False)
 for step in range(a.steps+1):
  with torch.no_grad():
   P=m.feats(m.decode(t));Z=(wt.transpose(1,2)@P).flatten(0,1)
   if step%a.inner_steps==0:
    left,sv,right=torch.linalg.svd((Z-dual)@metric_X.T,full_matrices=False);Y=(left@right)@X;target=Y+dual
   u,P,wt,Z,loss=fit(t,target)
  if step==a.steps:record(step,Y,Z,loss,False);break
  t.requires_grad_(True);uu=m.decode(t);ZZ=wt.transpose(1,2)@m.feats(uu);res=ZZ-target.reshape(R,C,-1);ll=(res*metric(res)).sum();ll.backward();grad=t.grad.detach();t=t.detach()
  with torch.no_grad():
   tf=t.flatten(0,1);v=torch.zeros((R*B,76),**m.d);v[:,0]=1;J=torch.zeros((R*B,75,76),**m.d)
   for j,(pp,cc) in enumerate(tree):
    co=(tf[:,j]*np.pi/128).cos();sn=(tf[:,j]*np.pi/128).sin();jp,jc=J[:,:,pp].clone(),J[:,:,cc].clone();J[:,:,pp]=co[:,None]*jp-sn[:,None]*jc;J[:,:,cc]=sn[:,None]*jp+co[:,None]*jc;vp,vc=v[:,pp].clone(),v[:,cc].clone();v[:,pp]=co*vp-sn*vc;v[:,cc]=sn*vp+co*vc;J[:,j,pp]=-np.pi/128*v[:,cc];J[:,j,cc]=np.pi/128*v[:,pp]
   h=(8*mean_metric*J.square().sum(-1).reshape(R,B,75)*wt.square().sum(-1)[...,None]).clamp(min=1e-18);order=(.5*h-grad.abs()).argsort(-1);direction=-grad.sign();candidates=[];ident=[];changes=[];residual=Z.reshape(R,C,-1)-target.reshape(R,C,-1);weighted=wt@metric(residual);wn=wt.square().sum(-1)
   global_words=[]
   if a.global_atoms:
    # With fixed coefficient row w_b and ||u||=1, the factor objective
    # differs from 2*u.T M_b*u by a constant. Its global minimizer is
    # the lowest eigenvector. Quantization and full refitting are exact.
    # L >= ||metric|| gives a quadratic upper model; ||u u.T||_F=1
    # makes its minimizer a lowest eigenvector even with a weighted metric.
    M=m.unpack(weighted-metric_lipschitz*wn[...,None]*P);_,vec=torch.linalg.eigh(M);direction_u=vec[...,0];direction_u*=torch.where((direction_u*u).sum(-1)>=0,1.,-1.)[...,None]
    for alpha in (.25,.5,1.):
     un=(1-alpha)*u+alpha*direction_u;un/=un.norm(dim=-1,keepdim=True);qw=encode(un.cpu().numpy().reshape(-1,76),tree,8,True);global_words.append(m.tensor(np.where(qw>=128,qw-256,qw).reshape(R,B,75)))
   for r in range(R):
    cp=[]
    for count in (1,2,4,8):
     ix=order[r,:,:count];delta=torch.zeros_like(t[r]);delta.scatter_(1,ix,direction[r].gather(1,ix));cp.append(t[r]+delta)
    cp.extend(q[r] for q in global_words)
    cp=torch.stack(cp,1);dp=m.feats(m.decode(cp))-P[r,:,None,:];pred=2*(dp*weighted[r,:,None,:]).sum(-1)+wn[r,:,None]*(dp*metric(dp)).sum(-1);pred[dp.square().sum(-1)<1e-22]=torch.inf;score,choice=pred.min(-1)
    for row in range(B):candidates.append(float(score[row]));ident.append((r,row));changes.append(cp[row,choice[row]]-t[r,row])
   ix=np.argsort(candidates);changed=False
   for size in (128,64,32,16,8,4,2,1):
    trial=t.clone()
    for j in ix[:size]:r,row=ident[j];trial[r,row]+=changes[j]
    tu,tp,tw,tz,tl=fit(trial,target)
    if tl<loss-1e-12:
     t,wt,Z,loss=trial,tw,tz,tl;accepted+=size;changed=True;break
   if (step+1)%a.inner_steps==0:dual+=a.dual_step*(Y-Z)
   if step%5==0:record(step,Y,Z,loss,changed)

if __name__=='__main__':main()
