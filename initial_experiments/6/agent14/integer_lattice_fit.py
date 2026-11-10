"""Integer row proposals from full projected angle Hessians.

The proposal relaxation uses all angle correlations. Only exact integer words
with an improving native signed rank-C refit are accepted.
"""
import sys,argparse,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
from lattice_babai import proposals as lattice_proposals
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);p.add_argument('--steps',type=int,default=100);p.add_argument('--sweeps',type=int,default=12);a=p.parse_args();torch.manual_seed(142130);z,u,w,tree,words=load(a.init);m=Metric(tree);torch.set_num_threads(2);R,B,C=w.shape;assert R==1;t=m.tensor(np.where(words>=128,words-256,words).reshape(R,B,75));wt=m.tensor(w);off=m.tensor(z['W']);start=time.time();history=[];stale=0;accepted=0
for step in range(a.steps):
 t.requires_grad_(True);u=m.decode(t);loss=m.loss(u,wt);loss.backward();g=t.grad[0].detach();t=t.detach();u=u.detach()
 with torch.no_grad():
  P=m.feats(u)[0];S=P@P.T;L=torch.linalg.cholesky(S);X=wt[0].T@P;D=X.T@X-m.G;Y=wt[0]@X;yn=Y.square().sum(-1);tf=t[0];v=torch.zeros((B,76),**m.d);v[:,0]=1;J=torch.zeros((B,75,76),**m.d)
  for j,(pp,cc) in enumerate(tree):
   co=(tf[:,j]*np.pi/128).cos();sn=(tf[:,j]*np.pi/128).sin();jp,jc=J[:,:,pp].clone(),J[:,:,cc].clone();J[:,:,pp]=co[:,None]*jp-sn[:,None]*jc;J[:,:,cc]=sn[:,None]*jp+co[:,None]*jc;vp,vc=v[:,pp].clone(),v[:,cc].clone();v[:,pp]=co*vp-sn*vc;v[:,cc]=sn*vp+co*vc;J[:,j,pp]=-np.pi/128*v[:,cc];J[:,j,cc]=np.pi/128*v[:,pp]
  gram=2*torch.bmm(J,J.transpose(1,2));H=gram.clone()
  for lo in range(0,B,64):
   hi=min(B,lo+64);cross=2*(J[lo:hi]@u[0].T)*(u[0,lo:hi]@u[0].T)[:,None,:];flat=cross.flatten(0,1);sol=torch.cholesky_solve(flat.T,L).T.reshape(hi-lo,75,B);H[lo:hi]-=cross@sol.transpose(1,2)
  H=(H+.0001*gram)*4*yn[:,None,None];diag=H.diagonal(dim1=1,dim2=2).clone();floor=diag.max(-1).values[:,None]*1e-7+1e-16;H.diagonal(dim1=1,dim2=2).add_(floor);diag=H.diagonal(dim1=1,dim2=2);proposals=[]
  relevance=(g.square()/diag.clamp(min=1e-30)).sum(-1);ids=relevance.argsort(descending=True)[:16].tolist();candidate_banks=[tf.clone() for _ in range(12)]
  for row in ids:
   cand=lattice_proposals(H[row].cpu().numpy(),g[row].cpu().numpy(),max_step=64)
   for j,(_,move) in enumerate(cand[:12]):candidate_banks[j][row]+=m.tensor(move)
  proposals=candidate_banks
  cp=torch.stack(proposals,1);dp=m.feats(m.decode(cp))-P[:,None,:];flat=dp.flatten(0,1);cross=flat@P.T;perp=(flat.square().sum(-1)-(cross*torch.cholesky_solve(cross.T,L).T).sum(-1)).clamp(min=0).reshape(B,-1);pred=4*(dp*(Y@D)[:,None,:]).sum(-1)+2*yn[:,None]*perp;pred[dp.square().sum(-1)<1e-22]=torch.inf;best,ki=pred.min(-1);chosen=cp[torch.arange(B,device=g.device),ki];order=best.argsort().tolist();err=float(loss.sqrt());changed=False;newerr=err
  groups=[order[:min(k,B)] for k in [B,128,64,32,16,8,4,2,1]]+[order[j:j+1] for j in range(1,10)]
  for group in groups:
   trial=t.clone();trial[0,group]=chosen[group];uu=m.decode(trial);ww,ne=m.project(uu,wt,1)
   if ne<err-1e-10:t=trial;wt=ww;newerr=ne;accepted+=len(group);changed=True;break
  stale=0 if changed else stale+1;row=dict(step=step,error=err,new_error=newerr,changed=changed,accepted_rows=accepted,seconds=time.time()-start);history.append(row);print(json.dumps(row),flush=True);save(f'parameters/{a.tag}.npz',m.decode(t),wt,tree,t,W=off.cpu().numpy());Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',arguments=vars(a),rows=history),indent=2)+'\n')
  if stale>=5:break
with torch.no_grad():
 m.project_gap(m.decode(t),wt,off);gg,ll,ss=m.physical(m.decode(t),wt,off);row=dict(final=float(m.loss(m.decode(t),wt).sqrt()),Lambda=float(ll),gap=float(gg),lambda_eff=float(ss.sqrt()),seconds=time.time()-start);print(json.dumps(row),flush=True);save(f'parameters/{a.tag}.npz',m.decode(t),wt,tree,t,W=off.cpu().numpy());Path(f'results/{a.tag}_final.json').write_text(json.dumps(row,indent=2)+'\n')
