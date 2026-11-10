"""Global integer proposals coupling angles on different dictionary rows."""
import sys,argparse,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
from lattice_qubo import solve
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);p.add_argument('--steps',type=int,default=100);p.add_argument('--variables',type=int,default=96);a=p.parse_args();torch.manual_seed(142134);z,u,w,tree,words=load(a.init);m=Metric(tree);torch.set_num_threads(2);R,B,C=w.shape;assert R==1;t=m.tensor(np.where(words>=128,words-256,words).reshape(R,B,75));wt=m.tensor(w);off=m.tensor(z['W']);start=time.time();rows=[];stale=0;accepted=0
for step in range(a.steps):
 t.requires_grad_(True);u=m.decode(t);loss=m.loss(u,wt);loss.backward();g=t.grad[0].detach();t=t.detach();u=u.detach()
 with torch.no_grad():
  P=m.feats(u)[0];S=P@P.T;L=torch.linalg.cholesky(S);X=wt[0].T@P;Y=wt[0]@X;YY=Y@Y.T;tf=t[0];v=torch.zeros((B,76),**m.d);v[:,0]=1;J=torch.zeros((B,75,76),**m.d)
  for j,(pp,cc) in enumerate(tree):
   co=(tf[:,j]*np.pi/128).cos();sn=(tf[:,j]*np.pi/128).sin();jp,jc=J[:,:,pp].clone(),J[:,:,cc].clone();J[:,:,pp]=co[:,None]*jp-sn[:,None]*jc;J[:,:,cc]=sn[:,None]*jp+co[:,None]*jc;vp,vc=v[:,pp].clone(),v[:,cc].clone();v[:,pp]=co*vp-sn*vc;v[:,cc]=sn*vp+co*vc;J[:,j,pp]=-np.pi/128*v[:,cc];J[:,j,cc]=np.pi/128*v[:,pp]
  raw=2*J.square().sum(-1);cross=2*(J@u[0].T)*(u[0]@u[0].T)[:,None,:];flat=cross.flatten(0,1);sol=torch.cholesky_solve(flat.T,L).T.reshape(B,75,B);diag=4*(raw-(cross*sol).sum(-1)+.0001*raw).clamp(min=1e-30)*YY.diag()[:,None];score=g.square()/diag.clamp(min=1e-30);ranking=score.flatten().argsort(descending=True).tolist();rrcounts={};ids=[]
  for i in ranking:
   r=i//75
   if rrcounts.get(r,0)>=2:continue
   ids.append(i);rrcounts[r]=rrcounts.get(r,0)+1
   if len(ids)==a.variables:break
  idx=torch.tensor(ids,device=P.device);rs=idx//75;ds=J.flatten(0,1)[idx];uv=u[0,rs];ju=ds@uv.T;rawgram=2*(ds@ds.T)*(uv@uv.T)+2*ju*ju.T;cr=flat[idx];so=sol.flatten(0,1)[idx];gram=rawgram-cr@so.T;H=4*(gram+.0001*rawgram)*YY[rs[:,None],rs[None,:]];gg=g.flatten()[idx];dd=torch.tensor([1.,-1.],**m.d).repeat(a.variables);HH=H.repeat_interleave(2,0).repeat_interleave(2,1)*dd[:,None]*dd[None,:];gg2=gg.repeat_interleave(2)*dd;hn=HH.cpu().numpy();bn=(gg2+.5*HH.diag()).cpu().numpy();penalty=max(float(abs(hn).max()),1e-14)*100
  for j in range(a.variables):hn[2*j,2*j+1]+=penalty;hn[2*j+1,2*j]+=penalty
  err=float(loss.sqrt());best=err;beststate=None;preds=[]
  for trial in range(8):
   selected,energy=solve(hn,bn,restarts=48,cap=128,seed=142134+1000*step+trial);preds.append(energy)
   if energy>=-1e-14:continue
   delta=torch.zeros(75*B,**m.d)
   for q in selected:delta[ids[q//2]]+=1 if q%2==0 else -1
   if not delta.any():continue
   proposal=t+delta.reshape(1,B,75);uu=m.decode(proposal);ww,ne=m.project(uu,wt,1)
   if ne<best-1e-10:best=ne;beststate=(proposal,ww,int(delta.count_nonzero()))
  changed=beststate is not None
  if changed:t,wt,count=beststate;accepted+=count;stale=0
  else:stale+=1
  row=dict(step=step,error=err,new_error=best,changed=changed,accepted_angles=accepted,min_predicted=min(preds),seconds=time.time()-start);rows.append(row);print(json.dumps(row),flush=True);save(f'parameters/{a.tag}.npz',m.decode(t),wt,tree,t,W=off.cpu().numpy());Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',arguments=vars(a),rows=rows),indent=2)+'\n')
  if stale>=8:break
with torch.no_grad():
 m.project_gap(m.decode(t),wt,off);gg,ll,ss=m.physical(m.decode(t),wt,off);row=dict(final=float(m.loss(m.decode(t),wt).sqrt()),Lambda=float(ll),gap=float(gg),lambda_eff=float(ss.sqrt()),seconds=time.time()-start);print(json.dumps(row),flush=True);save(f'parameters/{a.tag}.npz',m.decode(t),wt,tree,t,W=off.cpu().numpy());Path(f'results/{a.tag}_final.json').write_text(json.dumps(row,indent=2)+'\n')
