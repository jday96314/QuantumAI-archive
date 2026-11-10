"""Classical prototype: Euclid with high remainders and two growing cofactors.

For x^N and S, consecutive remainders Rprev,R and cofactors Bprev,B obey
Rprev*B + R*Bprev = x^N. After k degree-one quotients, only R[k:N-k]
needs storage. The omitted low coefficients follow from this identity.
Two cofactor polynomials and the retained remainder occupy N+1 field words.
This prototype rejects Euclidean degree gaps; it is not a decoder circuit.
"""
import json,random
from functools import lru_cache
from pathlib import Path
from import_published import load
m=load();mul=m.gmul;inv=m.ginv

def product(a,b):
 if not a or not b:return []
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]^=mul(x,y)
 return c

def xor(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]^=x
 for i,x in enumerate(b):c[i]^=x
 return c

def divide(a,b):
 a=list(a)
 while a and not a[-1]:a.pop()
 q=[0]*max(0,len(a)-len(b)+1)
 while len(a)>=len(b):
  d=len(a)-len(b);v=mul(a[-1],inv(b[-1]));q[d]=v
  for j,x in enumerate(b):a[d+j]^=mul(v,x)
  while a and not a[-1]:a.pop()
 return q,a

def step(N,k,Bprev,B,C):
 if not C[-1]:raise ValueError('degree gap or termination')
 top=inv(B[-1]);sub=mul(B[-2],mul(top,top)) if k else 0
 a=mul(top,inv(C[-1]));b=mul(sub^mul(a,C[-2]),inv(C[-1]));q=[b,a]
 Bnew=xor(Bprev,product(q,B))
 assert Bnew[-1]==inv(C[-1])
 assert Bnew[-2]==mul(C[-2],mul(inv(C[-1]),inv(C[-1])))
 # Compute only the required high band, by upper triangular division.
 out=[0]*(len(C)-2)
 for l in reversed(range(len(out))):
  degree=2*k+1+l;v=0
  for i,co in enumerate(Bnew):
   j=degree-i-k
   if 0<=j<len(C):v^=mul(co,C[j])
  for h in range(1,min(k,len(out)-1-l)+1):v^=mul(B[k-h],out[l+h])
  out[l]=mul(v,top)
 return list(B),Bnew,out,q

def backward(N,k,Bprev,Bnew,Cnew):
 # k denotes the old cofactor degree. Recover the two consumed coefficients.
 a=mul(Bnew[-1],inv(Bprev[-1]));b=mul(Bnew[-2]^mul(a,Bprev[-2] if k else 0),inv(Bprev[-1]));q=[b,a]
 oldprev=xor(Bnew,product(q,Bprev))[:k]
 lead=inv(Bnew[-1]);sub=mul(Bnew[-2],mul(lead,lead))
 old=[0]*(len(Cnew)+2);old[-2:]=[sub,lead]
 for l in reversed(range(len(Cnew))):
  degree=2*k+1+l;v=0
  for i,co in enumerate(Bprev):
   j=degree-i-(k+1)
   if 0<=j<len(Cnew):v^=mul(co,Cnew[j])
  for h in range(1,min(k+1,len(old)-1-l)+1):v^=mul(Bnew[k+1-h],old[l+h])
  old[l]=mul(v,inv(Bnew[-1]))
 return oldprev,list(Bprev),old

@lru_cache(None)
def triangular_products(n,degree):
 if n<=1 or degree<=0:return 0
 lo=n//2;hi=n-lo;l=min(lo,degree);h=min(hi,degree)
 pairs=sum(l+i-j<=degree for i in range(h) for j in range(l))
 return triangular_products(lo,degree)+triangular_products(hi,degree)+min(pairs,l+h-1)

def resource_lower_bound(N):
 rows=[]
 for k in range(N//2):
  n=N-2*k;remain=n-2
  off=triangular_products(remain,k+1)+triangular_products(remain,k)
  scales=remain+(remain if k else 0)
  boundary=sum(1 for l in range(remain) for h in [remain-l,remain+1-l] if h<=k+1)
  cofactor=min(max(0,2*k-1),k+2)
  rows.append(dict(k=k,remainder_words=n,off_diagonal_products=off,scalings=scales,boundary_products=boundary,cofactor_products=cofactor))
 total=sum(sum(r[z] for z in ['off_diagonal_products','scalings','boundary_products','cofactor_products']) for r in rows)
 return dict(scope='field-product estimate only; excludes pivot inverses, quotient conversion, controls, and phase kernel',N=N,rows=rows,field_products=total,T_field_products=42*total,storage_field_words=N+1)

def test():
 rng=random.Random(817818);cases=[];rejected=0
 for N,trials in [(4,100),(6,100),(10,100),(20,30),(36,10),(70,10)]:
  good=0
  while good<trials:
   S=[rng.randrange(4096) for _ in range(N)];Bp=[];B=[1];C=list(S);rp=[0]*N+[1];r=list(S);fp=[];f=[1];states=[]
   try:
    for k in range(N//2):
     if len(r)!=N-k or not r[-1]:raise ValueError('degree gap')
     states.append((Bp[:],B[:],C[:]))
     qb,rb=divide(rp,r)
     if len(qb)!=2:raise ValueError('degree gap')
     Bp,B,C,q=step(N,k,Bp,B,C)
     assert q==qb
     fp,f=f,xor(fp,product(qb,f));rp,r=r,rb
     assert B==f and Bp==fp
     assert C==(r+[0]*N)[k+1:N-k-1]
     old=backward(N,k,Bp,B,C);assert old==states[-1]
    for k in reversed(range(N//2)):
     Bp,B,C=backward(N,k,Bp,B,C);assert (Bp,B,C)==states[k]
    assert Bp==[] and B==[1] and C==S
   except ValueError:rejected+=1;continue
   good+=1
  cases.append(dict(N=N,trials=trials))
 result=dict(status='PASS',scope='generic degree-one Euclid path; exact forward, inverse, full-Euclid comparison',cases=cases,rejected_degree_gaps=rejected)
 Path('results/compact_euclid_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
 profile=resource_lower_bound(70);Path('results/compact_euclid_resources.json').write_text(json.dumps(profile,indent=2)+'\n');print(json.dumps({k:v for k,v in profile.items() if k!='rows'}))

if __name__=='__main__':test()
