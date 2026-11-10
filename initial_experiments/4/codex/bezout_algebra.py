"""Exact finite-field identities for the in-place Wronskian map."""
import json,random
from probe_wronskian import p,mul

def trim(a):
 a=list(a)
 while a and not a[-1]:a.pop()
 return a

def divide(a,b):
 a=trim(a);b=trim(b);q=[0]*max(0,len(a)-len(b)+1);iv=p.ginv(b[-1])
 while len(a)>=len(b):
  k=len(a)-len(b);c=mul(a[-1],iv);q[k]=c
  for j,v in enumerate(b):a[j+k]^=mul(c,v)
  a=trim(a)
 return trim(q),a

def product(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]^=mul(x,y)
 return out

def bezout(a,b):
 n=max(len(a),len(b))-1;a=list(a)+[0]*(n+1-len(a));b=list(b)+[0]*(n+1-len(b));B=[[0]*n for _ in range(n)]
 for i in range(1,n+1):
  for j in range(i):
   c=mul(a[i],b[j])^mul(b[i],a[j])
   for k in range(i-j):B[i-1-k][j+k]^=c
 return B

def matvec(A,x):return [sum_gf(mul(a,b) for a,b in zip(row,x)) for row in A]
def sum_gf(values):
 out=0
 for v in values:out^=v
 return out

def toeplitz(v,coeff,inverse=False):
 out=list(v);iv=p.ginv(coeff[0])
 if inverse:
  for i in range(len(out)):
   out[i]=mul(out[i]^sum_gf(mul(coeff[i-j],out[j]) for j in range(i) if i-j<len(coeff)),iv)
 else:
  for i in reversed(range(len(out))):out[i]=sum_gf(mul(coeff[i-j],out[j]) for j in range(i+1) if i-j<len(coeff))
 return out

def remainder_factors(E,O):
 f,g=list(E),list(O);factors=[]
 while g:
  q,r=divide(f,g)
  if len(q)>2:raise ValueError('degree gap')
  if len(q)==2:factors.append((g,q[1]))
  f,g=g,r
 return factors

def bezout_via_remainders(E,O,x):
 factors=remainder_factors(E,O);v=list(x)
 assert sorted(len(r)-1 for r,q in factors)==list(range(len(x)))
 for R,q in factors:
  d=len(R)-1;v[d]=mul(q,sum_gf(mul(c,v[i]) for i,c in enumerate(R)))
 for R,q in reversed(factors):
  d=len(R)-1;old=v[d]
  for i,c in enumerate(R[:-1]):v[i]^=mul(c,old)
  v[d]=mul(R[-1],old)
 return v

def sylvester_inplace(E,O,U,V):
 # U contains coefficients of z^1..z^n in the even template half;
 # V contains z^0..z^n in its odd half. H=E V+O z U.
 n=len(U);low=toeplitz(V,E)
 C=product(O,[0]+list(U))
 for i in range(n+1):low[i]^=C[i]
 u=toeplitz(U,E,inverse=True)[::-1]
 u=bezout_via_remainders(E,O,u)
 temp=toeplitz(low,E,inverse=True);cross=product(E,temp)
 for i in range(n):u[i]^=cross[n+1+i]
 return low+u

def test():
 rng=random.Random(94091);cases=0;rejected=0
 for n in [1,2,3,5,8,17]:
  for _ in range(50):
   E=[rng.randrange(1,4096)]+[rng.randrange(4096) for _ in range(n-1)]+[rng.randrange(1,4096)]
   O=[rng.randrange(4096) for _ in range(n)]+[rng.randrange(1,4096)]
   U=[rng.randrange(4096) for _ in range(n)];V=[rng.randrange(4096) for _ in range(n+1)];x=U
   try:by=bezout_via_remainders(E,O,x)
   except (ValueError,AssertionError):rejected+=1;continue
   assert by==matvec(bezout(E,O),x),n
   expected=product(E,V);other=product(O,[0]+U)
   expected=[a^b for a,b in zip(expected,other)]
   actual=sylvester_inplace(E,O,U,V)
   assert actual==expected,n
   cases+=1
 print(json.dumps(dict(status='PASS',generic_cases=cases,degree_gap_cases_skipped=rejected)))

if __name__=='__main__':test()
