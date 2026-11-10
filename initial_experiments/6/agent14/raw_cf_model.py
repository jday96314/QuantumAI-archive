"""Raw quadratic-node continued fraction with a projective infinity jet.

No scalar normalization or node swapping is required in the head extraction.
The two leading terms at infinity determine a finite affine terminal tail;
a completed infinite terminal tail is represented projectively.
"""
import json,random,time
from pathlib import Path
from quadratic_schur_model import F,node_maps,transform,bmul,einv
from derivative_inverse_model import trim,divide


def extract(s):
 n=len(s)//2;nodes,kap,P,maps=node_maps(n);v=transform(s,maps);heads=[];ghost=[]
 for i in range(n):
  b=v[i]>>12;a=(v[i]&4095)^F.mul(nodes[i]&4095,b);heads.append((a,b))
  for j in range(i+1,n):
   den=v[j]^a^bmul(nodes[j],b);ghost.append(den);v[j]=bmul(einv(den),kap[i]^kap[j]) if den else 0
 return heads,ghost


def endpoint(heads):
 p0=p1=q1=0;q0=1
 for a,b in heads:
  p0,p1,q0,q1=q0,q0^q1,p0^F.mul(b,q0),p1^F.mul(b,q1)^F.mul(a,q0)
 if q0:
  qi=F.inv(q0);h=F.mul(p0,qi);k=F.mul(p1,qi)^F.mul(h,F.mul(q1,qi));return [k,h],[1],None,(p0,p1,q0,q1)
 if not q1:return [1],[],None,(p0,p1,q0,q1)
 return [],[],'infinite_simple_pole',(p0,p1,q0,q1)


def synthesize(heads):
 N,B,exc,jet=endpoint(heads)
 if exc:return N,B,exc,jet
 _,kap,_,_=node_maps(len(heads))
 for i,(a,b) in reversed(list(enumerate(heads))):N,B=trim(F.add(F.product([a,b],N),F.product([kap[i],1,1],B))),N
 if not B:return N,B,'zero_denominator',jet
 if len(N)>=len(B):return N,B,'improper_final',jet
 return N,B,None,jet


def test(N=70,cases=4096):
 n=N//2;_,kap,P,_=node_maps(n);rng=random.Random(142111+N);rows={};bad=[];start=time.time();forced=0
 for case in range(cases):
  d=case%(n+1);pos=rng.sample(range(1,4096),d);val=[rng.randrange(1,4096) for _ in pos]
  if case%2 and d>=2:
   coeff=[F.mul(F.mul(x,F.evaluate(P,x)),F.inv(F.mul(x,x)^x^kap[0])) for x in pos];last=0
   for c,e in zip(coeff[:-1],val[:-1]):last^=F.mul(c,e)
   val[-1]=F.mul(last,F.inv(coeff[-1]));
   if not val[-1]:continue
   forced+=1
  heads,_=extract(F.syndrome(pos,val,N));A,B,exc,jet=synthesize(heads)
  if exc:rows[exc]=rows.get(exc,0)+1;continue
  label=f'degree_{len(B)-1}';rows[label]=rows.get(label,0)+1
  for x,e in zip(pos,val):
   if F.evaluate(B,x) or F.mul(F.evaluate(A,x),F.inv(F.evaluate(F.derivative(B),x)))!=F.mul(e,F.mul(x,F.evaluate(P,x))):bad.append(dict(case=case,d=d,root=x));break
  # Padding is squarefree iff B has coprime derivative, including irreducible
  # quadratic padding roots outside the scanned base field.
  a,b=B,F.derivative(B)
  while b:_,r=divide(a,b);a,b=b,r
  if len(a)>1:rows['nonsquarefree']=rows.get('nonsquarefree',0)+1
 r=dict(status='MODEL_PASS' if not bad else 'FAIL',N=N,cases=cases,forced=forced,rows=rows,bad=bad[:8],seconds=time.time()-start);Path(f'results/raw_cf_n{N}_model.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);assert not bad
if __name__=='__main__':
 for N,cases in [(6,2048),(14,2048),(70,4096)]:test(N,cases)
