"""Fixed-size integer dictionary exchanges, with exact PSD rank-C refits."""
import sys,argparse,json,time
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit','initial_experiments/5/chem-census']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
from census_common import initial
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);p.add_argument('--rounds',type=int,default=15);p.add_argument('--batch',type=int,default=8);a=p.parse_args();torch.manual_seed(142135);z,un,wn,tree,words=load(a.init);m=Metric(tree);torch.set_num_threads(2);u=m.tensor(un);w=m.tensor(wn);off=m.tensor(z['W']);R,B,C=w.shape;assert R==1;_,_,fu,fw,fo,fwords=initial();ix=np.load('../../5/chem-census/parameters/projected_B671.npz')['indices'];missing=np.setdiff1d(np.arange(719),ix);start=time.time();rows=[]
for iteration in range(a.rounds):
 with torch.no_grad():
  P=m.feats(u)[0];X=w[0].T@P;D=m.G-X.T@X;old_error=float(D.norm());Q,_=torch.linalg.qr(P.T);omega=torch.randn((len(D),48),**m.d);V,_=torch.linalg.qr(D@omega)
  for _ in range(3):V,_=torch.linalg.qr(D@V)
  ee,vv=torch.linalg.eigh(V.T@D@V);modes=(V@vv[:,-12:]).T;em,vm=torch.linalg.eigh(m.unpack(modes));ii=em.abs().argsort(-1,descending=True)[:,:8];cand=vm.gather(2,ii[:,None,:].expand(-1,76,-1)).transpose(1,2).flatten(0,1);cw=np.unique(encode(cand.cpu().numpy(),tree,8,True),axis=0);cw=np.concatenate([fwords[missing],cw]);cu=m.decode(m.tensor(cw));PP=m.feats(cu);perp=PP-(PP@Q)@Q.T;norm=perp.square().sum(-1);v=perp/norm.sqrt().clamp(min=1e-12)[:,None];dv=v@D;gain=2*(dv@Q).square().sum(-1)+(dv*v).sum(-1).clamp(min=0).square();gain[norm<1e-12]=-torch.inf;existing={tuple(q) for q in words}
  for j,row in enumerate(cw):
   if tuple(row) in existing:gain[j]=-torch.inf
  chosen=[];basis=[]
  for j in gain.argsort(descending=True).tolist():
   if not torch.isfinite(gain[j]):continue
   vj=perp[j].clone()
   if basis:
    qb=torch.stack(basis);vj-=(qb@vj)@qb
   if float(vj.norm())<1e-6:continue
   chosen.append(j);basis.append(vj/vj.norm())
   if len(chosen)>=a.batch:break
  if not chosen:break
  uu=torch.cat([u,cu[chosen][None]],1);ww=torch.nn.functional.pad(w,(0,0,0,len(chosen)));ww,grown_error=m.project(uu,ww,1);aw=np.concatenate([words,cw[chosen]]);deleted=[]
  for _ in chosen:
   ps=m.feats(uu)[0];S=ps@ps.T;si=torch.linalg.inv(S);K=ww[0]@ww[0].T;d=si.diag();importance=2*(K@S@K).diag()/d-(K.diag()/d).square();best=1e10;keepbest=None
   for remove in importance.argsort()[:6].tolist():
    keep=[j for j in range(len(ps)) if j!=remove];pu=uu[:,keep];pw,er=m.project(pu,ww[:,keep],1)
    if er<best:best=er;keepbest=(pu,pw,aw[keep],remove)
   uu,ww,aw,remove=keepbest;deleted.append(remove)
  new_error=float(m.loss(uu,ww).sqrt());accepted=new_error<old_error-1e-10
  if accepted:u,w,words=uu,ww,aw
  row=dict(iteration=iteration,old_error=old_error,grown_error=grown_error,exchange_error=new_error,accepted=accepted,added=len(chosen),seconds=time.time()-start);rows.append(row);print(json.dumps(row),flush=True);save(f'parameters/{a.tag}.npz',u,w,tree,words,W=off.cpu().numpy());Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',arguments=vars(a),rows=rows),indent=2)+'\n')
with torch.no_grad():
 m.project_gap(u,w,off);g,l,s=m.physical(u,w,off);row=dict(final=float(m.loss(u,w).sqrt()),Lambda=float(l),gap=float(g),lambda_eff=float(s.sqrt()),seconds=time.time()-start);print(json.dumps(row),flush=True);save(f'parameters/{a.tag}.npz',u,w,tree,words,W=off.cpu().numpy());Path(f'results/{a.tag}_final.json').write_text(json.dumps(row,indent=2)+'\n')
