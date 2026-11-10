"""Integer geometry training against the exact optimal rank-C PSD fit.

A straight-through derivative is only a proposal heuristic. The forward
objective and all saved checkpoints are exact integer-angle tensor fits.
"""
import sys,argparse,json,time,math
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
from robust_fixed_fit import continuous_encode
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);p.add_argument('--steps',type=int,default=6000);p.add_argument('--lr',type=float,default=.01);p.add_argument('--checkpoint',type=int,default=200);a=p.parse_args();torch.manual_seed(142133);z,u,w,tree,words=load(a.init);m=Metric(tree);torch.set_num_threads(2);R,B,C=w.shape;assert R==1;warm=m.tensor(w);off=m.tensor(z['W']);coord=m.tensor(continuous_encode(u,tree)).requires_grad_();eye=torch.eye(B,**m.d);opt=torch.optim.Adam([coord],lr=a.lr);const=m.G.square().sum();start=time.time();rows=[];best=1e10

def state():
 word=coord+(coord.round()-coord).detach();un=m.decode(word);P=m.feats(un)[0];L=torch.linalg.cholesky(P@P.T+1e-12*eye);H=torch.linalg.solve_triangular(L,P@m.G@P.T,upper=False);H=torch.linalg.solve_triangular(L,H.T,upper=False).T;ev=torch.linalg.eigvalsh((H+H.T)/2);err=(const-ev[-C:].clamp(min=0).square().sum()).clamp(min=1e-20);return un,err
for step in range(a.steps+1):
 if step%a.checkpoint==0 or step==a.steps:
  with torch.no_grad():
   un,err=state();wp,pe=m.project(un,warm,1);assert abs(pe**2-float(err))<1e-7,(pe,float(err));oo=off.clone();m.project_gap(un,wp,oo);g,l,s=m.physical(un,wp,oo);row=dict(step=step,integer_error=pe,Lambda=float(l),gap=float(g),lambda_eff=float(s.sqrt()),weight_norm=float(wp.norm()),seconds=time.time()-start);rows.append(row);print(json.dumps(row),flush=True)
   if pe<best:best=pe;save(f'parameters/{a.tag}_best.npz',un,wp,tree,coord.round(),W=oo.cpu().numpy())
   Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',method='exact rank-C PSD optimum at integer geometry with STE proposal gradients',arguments=vars(a),rows=rows,finished=step==a.steps),indent=2)+'\n')
 if step==a.steps:break
 rate=.5*(1+math.cos(math.pi*step/a.steps))*.99+.01;opt.param_groups[0]['lr']=a.lr*rate;opt.zero_grad();un,err=state();err.backward();torch.nn.utils.clip_grad_norm_([coord],1.);opt.step()
