"""Exact blocked C-fraction tail identity, a new architecture experiment.

Reference recurrence: https://dlmf.nist.gov/3.10 . All field and polynomial
operations use the existing independent round-five algebra, read-only.
"""
import sys,json,random
from pathlib import Path
sys.path.append(str(Path('../../5/dqi-core-phase').resolve()))
from algebra import Field
F=Field()

def fixed(p,n):return p[:n]+[0]*max(0,n-len(p))
def sprod(a,b,n):return fixed(F.product(a,b),n)
def quotient(a,b,n):return sprod(a,F.reciprocal(b,n),n)
def extract(s,count):
 c=[];s=list(s)
 for _ in range(count):
  a=s[0]
  if not a:raise ZeroDivisionError()
  c.append(a);s=[F.mul(a,v) for v in F.reciprocal(s,len(s))[1:]]
 return c,s

def transfer(c):
 A,B,C,D=[1],[],[],[1]
 for a in c:
  A,B,C,D=[0]+B,F.add(F.scale(A,a),B),[0]+D,F.add(F.scale(C,a),D)
 return A,B,C,D

def blocked(s,b):
 c,_=extract(s[:b],b);A,B,C,D=transfer(c);head=s[:b];t=s[b:];n=len(t)
 u=F.add(B,F.product(D,head))
 assert not any(u[:b])
 u=u[b:]
 N=fixed(F.add(sprod(D,t,n),u),n)
 K=1
 for a in c:K=F.mul(K,a)
 determinant=F.add(F.product(A,D),F.product(B,C))
 assert determinant==[0]*b+[K]
 L=C[1:];assert C[0]==0 and L[0]==D[0]==1
 # In-place arithmetic skeleton: unit multiplication, scale, conjugated
 # reciprocal with one virtual constant word, unit division, multiplication.
 v=F.scale(N,F.inv(K));v=sprod(L,v,n)
 v=F.reciprocal([1]+v,n+1)[1:]
 v=quotient(v,L,n);v=sprod(D,v,n)
 return c,v

def main():
 rng=random.Random(142026);rows=[]
 for n,b in [(14,2),(14,4),(22,6),(70,4),(70,8)]:
  cases=0;rejected=0
  while cases<32:
   s=[rng.randrange(4096) for _ in range(n)]
   try:c,w=extract(s,b);bc,v=blocked(s,b)
   except ZeroDivisionError:rejected+=1;continue
   assert c==bc and w==v,(n,b,c,bc,w,v)
   cases+=1
  rows.append(dict(N=n,block=b,cases=cases,zero_pivot_rejections=rejected))
 res=dict(status='ALGEBRA_PASS',scope='generic nonzero-pivot blocks only; no quantum resource or score claim',rows=rows)
 Path('results/blocked_cf_identity.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
if __name__=='__main__':main()
