"""Heuristic straight-through word training with exact integer checkpoints.

Every forward objective uses integer Givens words. The backward surrogate is
explicitly heuristic; saved models are judged by exact decoded native tensors.
"""
import sys,argparse,json,time,math
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
from robust_fixed_fit import continuous_encode
p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);p.add_argument('--steps',type=int,default=6000);p.add_argument('--lr',type=float,default=.0025);p.add_argument('--weight-lr',type=float,default=1e-5);p.add_argument('--lambda-reg',type=float,default=1e-7);p.add_argument('--checkpoint',type=int,default=200);a=p.parse_args();torch.manual_seed(142132);z,u,w,tree,words=load(a.init);m=Metric(tree);torch.set_num_threads(2);R,B,C=w.shape;ut=m.tensor(u);wt=m.tensor(w);off=m.tensor(z['W']);scale=40.;coord=m.tensor(continuous_encode(u,tree)/scale).requires_grad_();P=m.feats(ut);eye=torch.eye(B,**m.d);L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-10*eye);f=(L.transpose(1,2)@wt).detach().requires_grad_();opt=torch.optim.Adam([dict(params=[coord],lr=a.lr),dict(params=[f],lr=a.weight_lr)]);const=m.G.square().sum();start=time.time();rows=[];best=1e10

def state():
 master=coord*scale;word=master+(master.round()-master).detach();un=m.decode(word);P=m.feats(un);L=torch.linalg.cholesky(P@P.transpose(1,2)+1e-10*eye);wn=torch.linalg.solve_triangular(L.transpose(1,2),f,upper=True);X=(wn.transpose(1,2)@P).flatten(0,1);err=(const-2*(X*(X@m.G)).sum()+(X@X.T).square().sum()).clamp(min=1e-20);return un,wn,err,L
for step in range(a.steps+1):
 if step%a.checkpoint==0 or step==a.steps:
  with torch.no_grad():
   un,wn,err,L=state();actual=float(m.loss(un,wn).sqrt());wp,pe=m.project(un,wn,1);oo=off.clone();m.project_gap(un,wp,oo);g,l,s=m.physical(un,wp,oo);row=dict(step=step,integer_error=actual,refit_error=pe,Lambda=float(l),gap=float(g),lambda_eff=float(s.sqrt()),weight_norm=float(wp.norm()),seconds=time.time()-start);rows.append(row);print(json.dumps(row),flush=True)
   if pe<best:
    best=pe;save(f'parameters/{a.tag}_best.npz',un,wp,tree,(coord*scale).round(),W=oo.cpu().numpy())
   Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',method='straight-through backward surrogate, exact integer forward and native tensor refits',arguments=vars(a),rows=rows,finished=step==a.steps),indent=2)+'\n')
 if step==a.steps:break
 rate=.5*(1+math.cos(math.pi*step/a.steps))*.99+.01;opt.param_groups[0]['lr']=a.lr*rate;opt.param_groups[1]['lr']=a.weight_lr*rate;opt.zero_grad();un,wn,err,L=state();smooth=lambda x:torch.sqrt(x.square()+1e-14);h2=.25*(smooth(wn).sum(1)+smooth(off+wn.sum(1))).square().sum();loss=err+a.lambda_reg*h2;loss.backward();torch.nn.utils.clip_grad_norm_([coord],1.);opt.step()
