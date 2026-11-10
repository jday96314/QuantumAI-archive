"""Small resultant formulas for a head stored as denominator continuants.

A(0)=D(0)=1, deg A=b, deg D=b-1, determinant scalar p. The affine tail
P obeys A*Q + x^2*D*P = p, with deg Q<=b, deg P<b. Cramer's rule gives
P_i=p*N_i/resultant(A,D), avoiding reconstruction of all 2b moments.
Formulas are characteristic-two exact polynomials, not numerical guesses.
"""
from functools import lru_cache
from continuant_polynomials import add,mul,Compiler


@lru_cache(None)
def derive(b):
 n=2*b-1;z=frozenset();one=frozenset({(0,)*n});v=[frozenset({tuple(int(i==j) for j in range(n))}) for i in range(n)]
 A=[one,*v[:b]];D=[one,*v[b:]]
 def get(p,i):return p[i] if 0<=i<len(p) else z
 matrix=[];rhs=[]
 for k in range(2,2*b+1):
  matrix.append([*[get(A,k-j) for j in range(2,b+1)],*[get(D,k-2-i) for i in range(b)]])
  rhs.append(add(get(A,k),mul(A[1],get(A,k-1))))
 def det(M):
  @lru_cache(None)
  def rec(mask):
   row=mask.bit_count()
   if row==n:return one
   return add(*[mul(M[row][j],rec(mask|(1<<j))) for j in range(n) if not(mask>>j&1)])
  return rec(0)
 denominator=det(matrix);numerators=[]
 for col in range(b-1,n):
  M=[list(row) for row in matrix]
  for i in range(n):M[i][col]=rhs[i]
  numerators.append(det(M))
 return denominator,numerators


@lru_cache(None)
def compiled(b,precision=None):
 denominator,numerators=derive(b)
 if precision is not None and precision>b:
  count=2*(precision-b);n=2*b-1;z=frozenset();one=frozenset({(0,)*n})
  A=[one,*[frozenset({tuple(int(i==j) for j in range(n))}) for i in range(b)]]
  R=[]
  for k in range(count):R.append(add(numerators[k] if k<len(numerators) else z,*[mul(A[j],R[k-j]) for j in range(1,min(b,k)+1)]))
  numerators=[add(numerators[count+i] if count+i<len(numerators) else z,*[mul(A[j],R[count+i-j]) for j in range(b+1) if 0<=count+i-j<count]) for i in range(b)]
 c=Compiler(2*b-1);roots=[c.compile(denominator),*[c.compile(p) for p in numerators]];return c,roots


def value(m,p,vars):
 out=0
 for term in p:
  a=1
  for x,e in zip(vars,term):a=m.gmul(a,m.gpow(int(x),e))
  out^=a
 return out

if __name__=='__main__':
 import argparse,json,time
 from pathlib import Path
 import numpy as np
 from import_published import load
 import check_general_jacobi as ref
 p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,nargs='+',default=[2,3,4]);args=p.parse_args();m=load();rng=np.random.default_rng(2026102604);rows=[]
 for b in args.pairs:
  start=time.time();den,num=derive(b);cases=0
  for _ in range(256):
   H=rng.integers(1,4096,2*b).tolist();A,D,p,P=ref.parameters(m,H);vs=[*A[1:],*D[1:]];de=value(m,den,vs);assert de
   expect=[m.gmul(m.gmul(value(m,x,vs),p),m.ginv(de)) for x in num];assert P==expect,(b,P,expect);cases+=1
  c,roots=compiled(b);ops=set().union(*(c.operations[r] for r in roots));row=dict(pairs=b,cases=cases,denominator_monomials=len(den),numerator_monomials=[len(x) for x in num],shared_field_products=len(ops),nodes=len(c.nodes),seconds=time.time()-start,status='PASS_SYMBOLIC_AFFINE_IDENTITY')
  rows.append(row);print(json.dumps(row),flush=True);Path('results/continuant_frame_affine_formulas.json').write_text(json.dumps(rows,indent=2)+'\n')
