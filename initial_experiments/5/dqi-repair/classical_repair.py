"""Independent exact native-field lookahead Euclid and six-moment gap repair.

A missing first pivot needs a degree-two quotient; replacing it with one is
not a valid quotient update. Two adjacent Padé coordinates recover the state
without division by that pivot. This file does not claim a global quantum adapter.
"""
import json
from pathlib import Path
import numpy as np
EX=np.fromfile('parameters/antilog.bin',dtype=np.int32).tolist();LG=np.fromfile('parameters/log.bin',dtype=np.int32).tolist()
def mul(a,b):return EX[(LG[a]+LG[b])%4095] if a and b else 0
def inv(a):return EX[-LG[a]%4095] if a else 0
def trim(a):
 a=list(a)
 while a and not a[-1]:a.pop()
 return a
def product(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]^=mul(x,y)
 return trim(c)
def plus(a,b):
 c=list(a)+[0]*max(0,len(b)-len(a))
 for i,v in enumerate(b):c[i]^=v
 return trim(c)
def divide(a,b):
 a=trim(a);b=trim(b);q=[0]*max(0,len(a)-len(b)+1)
 while len(a)>=len(b):
  d=len(a)-len(b);v=mul(a[-1],inv(b[-1]));q[d]=v
  for j,x in enumerate(b):a[d+j]^=mul(v,x)
  a=trim(a)
 return trim(q),a
def monic(a):return [mul(v,inv(a[-1])) for v in a]
def repair(F):
 e,d,a,b,c,f=F;D=mul(d,d)^mul(e,a)
 if not D:raise ValueError('singular second minor')
 u=mul(mul(d,b)^mul(a,a),inv(D));v=mul(mul(e,b)^mul(d,a),inv(D))
 p=c^mul(a,u)^mul(b,v)
 if not p:raise ValueError('singular third minor')
 al=mul(f^mul(d,mul(u,u))^mul(b,mul(v,v)),inv(p));z=mul(p,inv(D))
 return [u,v,1],[mul(al,u)^mul(z,d),u^mul(al,v)^mul(z,e),v^al,1],inv(p)
def undo(U,V,pi):
 u,v=U[:2];v0,v1,v2=V[:3];al=v2^v;r0=v0^mul(al,u);r1=v1^u^mul(al,v)
 z=mul(mul(r0,r0)^mul(mul(r1,r1),u)^mul(mul(r0,r1),v),pi)
 e=mul(r1,inv(z));d=mul(r0,inv(z));a=mul(e,u)^mul(d,v);b=mul(d,u)^mul(a,v)
 return [e,d,a,b,inv(pi)^mul(a,u)^mul(b,v),mul(al,inv(pi))^mul(d,mul(u,u))^mul(b,mul(v,v))]
def euclid(F):
 P=[0]*len(F)+[1];R=trim(list(reversed(F)));bp=[];b=[1];out={};pattern=[]
 while R:
  q,Rn=divide(P,R);bn=plus(bp,product(q,b));P,R=R,Rn;bp,b=b,bn;pattern.append(len(q)-1)
  out[len(b)-1]=monic(b)
  if len(b)>3:break
 return out,pattern
if __name__=='__main__':
 rng=np.random.default_rng(6841729);counts={};good=0;skip=0
 for i in range(20000):
  F=rng.integers(0,4096,6).tolist()
  if i%2:F[0]=0
  try:U,V,pi=repair(F)
  except ValueError:skip+=1;continue
  assert undo(U,V,pi)==F
  out,pat=euclid(F);assert out[2]==U and out[3]==V,(F,out,U,V)
  key=','.join(map(str,pat));counts[key]=counts.get(key,0)+1;good+=1
 result=dict(status='PASS',scope='classical exact local degree-gap correction, inverse and full Euclid comparison; no global reversible decoder claim',cases=good,singular_rejected=skip,degree_patterns=counts)
 Path('results/classical_local_repair.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
