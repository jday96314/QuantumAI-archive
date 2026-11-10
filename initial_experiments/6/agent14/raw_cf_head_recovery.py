"""Check that each consumed raw head is recovered by its node quotient."""
import json,random
from pathlib import Path
from raw_cf_model import extract,endpoint
from quadratic_schur_model import F,node_maps,bmul,einv,emul
from derivative_inverse_model import trim
rng=random.Random(142113);rows=[]
for N in [6,14,70]:
 n=N//2;nodes,kap,P,_=node_maps(n);fail=[];denzero=0;cases=512
 def ev(poly,x):
  out=0
  for c in reversed(poly):out=emul(out,x)^c
  return out
 for case in range(cases):
  d=case%(n+1);pos=rng.sample(range(1,4096),d);val=[rng.randrange(1,4096) for _ in pos];heads,_=extract(F.syndrome(pos,val,N));A,B,exc,jet=endpoint(heads)
  if exc:continue
  for i,(a,b) in reversed(list(enumerate(heads))):
   A,B=trim(F.add(F.product([a,b],A),F.product([kap[i],1,1],B))),A
   dv=ev(B,nodes[i]);nv=ev(A,nodes[i]);q=emul(nv,einv(dv)) if dv else 0;want=a^bmul(nodes[i],b)
   if not dv:denzero+=1
   if q!=want:fail.append((case,d,i))
 r=dict(N=N,cases=cases,zero_divisors=denzero,errors=fail[:12],error_count=len(fail));rows.append(r);print(json.dumps(r),flush=True)
Path('results/raw_cf_head_recovery.json').write_text(json.dumps(rows,indent=2)+'\n')
