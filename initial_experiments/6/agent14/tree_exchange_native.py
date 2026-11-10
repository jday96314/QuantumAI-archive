"""Exact native-word NNI tree search on the repaired B671 dictionary.

Read-only NNI generation and signed quantization from round-5 chem-tree;
all numerical tensor checks and outputs remain local.
"""
import sys,argparse,json,time,importlib.util
from pathlib import Path
ROOT=Path.cwd().parents[2]
for p in ['initial_experiments/5/chem-qfit','initial_experiments/5/chem-audit']:sys.path.append(str(ROOT/p))
sys.path.insert(0,str(Path.cwd()))
from audit_core import *
spec=importlib.util.spec_from_file_location('r5_tree_ops',ROOT/'initial_experiments/5/chem-tree/trees.py');ops=importlib.util.module_from_spec(spec);spec.loader.exec_module(ops)
p=argparse.ArgumentParser();p.add_argument('--init',default='parameters/lll_B671_native.npz');p.add_argument('--tag',default='tree_B671');p.add_argument('--rounds',type=int,default=3);p.add_argument('--support',type=int,default=8);a=p.parse_args();z,un,wn,tree,words=load(a.init);m=Metric(tree);torch.set_num_threads(2);u=m.tensor(un);w=m.tensor(wn);off=m.tensor(z['W']);rows=[];start=time.time();best=float(m.loss(u,w).sqrt());initial=best
for repeat in range(a.rounds):
 accepted=None;roundbest=best;proposals=list(ops.nni_neighbors(tree,a.support));print(json.dumps(dict(repeat=repeat,proposals=len(proposals),baseline=best)),flush=True)
 for index,tr in enumerate(proposals):
  nw,nu=ops.quantize(u.cpu().numpy(),tr);tu=m.tensor(nu);ww,error=m.project(tu,w,1);row=dict(repeat=repeat,index=index,error=error,baseline=best,seconds=time.time()-start);rows.append(row)
  if error<roundbest-1e-10:roundbest=error;accepted=(tr,nw,tu,ww);save(f'parameters/{a.tag}_best.npz',tu,ww,tr,nw,W=off.cpu().numpy());print(json.dumps(dict(**row,new_best=True)),flush=True)
  if index%20==0:print(json.dumps(row),flush=True)
  Path(f'results/{a.tag}.json').write_text(json.dumps(dict(status='NUMERICAL_SCREEN',arguments=vars(a),initial=initial,best=roundbest,rows=rows),indent=2)+'\n')
 if accepted is None:break
 tree,words,u,w=accepted;best=roundbest
with torch.no_grad():
 m.project_gap(u,w,off);g,l,s=m.physical(u,w,off);result=dict(status='NUMERICAL_SCREEN',initial=initial,final=best,Lambda=float(l),gap=float(g),lambda_eff=float(s.sqrt()),seconds=time.time()-start);save(f'parameters/{a.tag}.npz',u,w,tree,words,W=off.cpu().numpy());Path(f'results/{a.tag}_final.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
