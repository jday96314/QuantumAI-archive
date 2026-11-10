"""Three-field infinity jet with a tail gauge making the final B monic."""
import json,random
from pathlib import Path
from raw_cf_model import extract
from quadratic_schur_model import F,node_maps
from derivative_inverse_model import trim


def endpoint(heads):
 p=d=0;q=1
 for a,b in heads:p,q,d=q,p^F.mul(b,q),d^F.mul(p,q)^F.mul(a^b,F.mul(q,q))
 if q:return [F.mul(d,F.inv(q)),p],[q],None,(p,q,d)
 if not d:return [p],[],None,(p,q,d)
 return [],[],'infinite_simple_pole',(p,q,d)


def synthesize(heads):
 A,B,e,j=endpoint(heads)
 if e:return A,B,e,j
 _,kap,_,_=node_maps(len(heads))
 for i,(a,b) in reversed(list(enumerate(heads))):A,B=trim(F.add(F.product([a,b],A),F.product([kap[i],1,1],B))),A
 return A,B,None,j


def test(N=70,cases=2048):
 rng=random.Random(142136+N);n=N//2;_,_,P,_=node_maps(n);counts={};bad=[]
 for case in range(cases):
  d=case%(n+1);pos=rng.sample(range(1,4096),d);val=[rng.randrange(1,4096) for _ in pos];heads,_=extract(F.syndrome(pos,val,N));A,B,e,j=synthesize(heads)
  if e:counts[e]=counts.get(e,0)+1;continue
  key=f'degree{len(B)-1}';counts[key]=counts.get(key,0)+1
  if not B or B[-1]!=1:bad.append(dict(case=case,lead=B[-1] if B else None));continue
  for x,v in zip(pos,val):
   if F.evaluate(B,x) or F.mul(F.evaluate(A,x),F.inv(F.evaluate(F.derivative(B),x)))!=F.mul(v,F.mul(x,F.evaluate(P,x))):bad.append(dict(case=case,root=x));break
 r=dict(status='MODEL_PASS' if not bad else 'FAIL',N=N,cases=cases,counts=counts,bad=bad[:8]);Path(f'results/raw_three_jet_n{N}_model.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);assert not bad
if __name__=='__main__':
 for N in [6,14,70]:test(N)
