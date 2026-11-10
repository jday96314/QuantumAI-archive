"""Reproduce read-only exact-word checkpoints and preserve matched refits."""
import json,hashlib,time
from core import *
from dfthcsos.verify import load_inputs,evaluate_words
from dfthcsos.scoreinfra.dfthc_math import build_h1_prime,reconstruct_x

prior=ROOT/'initial_experiments/5/chem-qfit'
tags=['incumbent_control','discrete_control','legacy_exact_split_B479',
      'discrete_claude_B479','qd_b447_C175','hardfreeze_B447','anneal_B671',
      'pruned_B639','greedy_B711','constrained_B711']
rows=[];m=None
for tag in tags:
    path=prior/f'candidates/{tag}.npz';z,u,w,tree,words=load(path)
    if m is None:m=Metric(tree)
    t=m.tensor(words.reshape(*u.shape[:2],75));ut=m.decode(t);wt=m.tensor(w)
    err=float(m.loss(ut,wt).sqrt());before=err
    save(HERE/f'candidates/start_{tag}.npz',ut,wt,tree,t,**({'W':z['W']} if 'W' in z and z['W'].ndim==2 else {}))
    trace=[];st=time.time()
    for sweep in range(101):
        if sweep in [0,1,4,20,100]:
            trace.append(dict(sweep=sweep,error_2=err,seconds=time.time()-st))
            save(HERE/f'candidates/refit{ sweep }_{tag}.npz',ut,wt,tree,t,err2=err)
        if sweep<100:wt,err=m.project(ut,wt,1)
    row=dict(tag=tag,shape=list(w.shape),brot=8,source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),initial_error=before,refits=trace)
    rows.append(row);print(json.dumps(row),flush=True)
(HERE/'results/reproduction.json').write_text(json.dumps(rows,indent=2))
inp=load_inputs(INC/'submission.npz');base=dict(np.load(INC/'submission.npz'))
h0=build_h1_prime(inp.ref['h1_thc'],inp.S,inp.beta_1,inp.ref['electrons'],np.zeros_like(inp.W),reconstruct_x(base['U'],inp.w))
np.savez(HERE/'candidates/native_constants.npz',h0=h0,constant=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1)
met=evaluate_words(inp,base['angle_words'],base['tree_edges'],8)
(HERE/'results/incumbent_physical.json').write_text(json.dumps(met,indent=2))
print(json.dumps(met),flush=True)
