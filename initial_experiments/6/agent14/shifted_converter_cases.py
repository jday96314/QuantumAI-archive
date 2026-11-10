"""Identify zero-constant Euclidean gap charts in fixed saved regressions."""
import json
from pathlib import Path
import numpy as np
from raw_cf_model import extract
from raw_three_jet_model import synthesize
from derivative_inverse_model import divide,trim
from quadratic_schur_model import F
rows=[]
def shift(p):
 out=[0]*len(p)
 for j,c in enumerate(p):
  for i in range(j+1):
   if i&~j==0:out[i]^=c
 return trim(out)
for path in sorted(Path('results/raw_cf_seedless_failures').glob('*.npz')):
 z=np.load(path)
 for s,idx in zip(z['syndromes'],z['shot_indices']):
  h,_=extract(s.tolist());A,B,e,j=synthesize(h)
  if len(B)==35:A=[0]+A;B=[0]+B
  records=[]
  for translated in [False,True]:
   C=shift(B) if translated else B;f=[F.pow(x,2048) for x in C[::2]];g=[F.pow(x,2048) for x in C[1::2]];_,r=divide(f,g);f,g=g,r;expected=17;used=False;issues=[]
   while expected:
    if len(g)!=expected:
     issues.append(dict(expected=expected,size=len(g),f0=f[0] if f else None,used=used))
     if not used and len(g)==expected-1:g=[0]+g;used=True
     else:break
    if expected==1:break
    _,r=divide(f,g);f,g=g,r;expected-=1
   records.append(dict(shift=translated,constant=C[0],issues=issues))
  rows.append(dict(shot=int(idx),records=records))
Path('results/shifted_converter_cases.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
