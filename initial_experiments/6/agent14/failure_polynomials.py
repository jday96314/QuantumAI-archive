"""Reconstruct exact polynomial inputs of saved common regression cases."""
import json
from pathlib import Path
import numpy as np
from schur_extract_reciprocal import model,controls
from quadratic_schur_model import node_maps,F
from schur_fold import reference
from extension_inverse import setup
m,_,fw,bk,tr=setup();files=sorted(Path('results/counterfree_failures').glob('*.npz'));zs=[dict(np.load(p)) for p in files];ss=np.concatenate([z['syndromes'] for z in zs]);ids=np.concatenate([z['shot_indices'] for z in zs]);_,kap,_,_=node_maps(35);rows=[];AA=[];BB=[];valid=[]
for s,idx in zip(ss,ids):
 heads,g,flags,ghost,pivot=model(s.tolist());N=[];B=[]
 try:
  for i in reversed(range(35)):
   kk=kap[i]
   for flag,delta in controls(70,i):
    if flags[flag]:kk^=F.mul(delta,delta)^delta
   N,B=reference(tr(heads[2*i],fw),tr(heads[2*i+1],fw),N,B,tr(kk,fw))
 except AssertionError:
  rows.append(dict(shot=int(idx),status='UNSUPPORTED_SCHUR'));continue
 c=next((i for i,v in enumerate(B+[1]) if v),35);rows.append(dict(shot=int(idx),status='POLYNOMIAL',padding=c));AA.append(N);BB.append(B);valid.append(idx)
np.savez('parameters/counterfree_failure_polynomials.npz',A=np.array(AA),B=np.array(BB),shot_indices=np.array(valid));Path('results/counterfree_failure_polynomials.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows),flush=True)
