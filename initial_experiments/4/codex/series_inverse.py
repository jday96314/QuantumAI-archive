"""Ancilla-light normalized formal-series inversion in characteristic two.

F^(2^k)=1 mod z^n once 2^k>=n. Compose powering by 3, 5, 17, ...
(the products are 2^(2^r)-1) to obtain F^-1. Each power map splits
F=A+z^m B and uses only triangular convolutions and a middle product.
The constant coefficient is implicitly one and must be one on input.
"""
from functools import lru_cache
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
ROOT=Path(__file__).resolve().parent

class SeriesArithmetic:
 def __init__(self,module,emitter):self.m=module;self.E=emitter
 @lru_cache(None)
 def rows(self,c,q=1):
  m=self.m;cols=[m.gmul(c,m.gpow(1<<i,q)) for i in range(12)]
  return tuple(sum(((v>>j)&1)<<i for i,v in enumerate(cols)) for j in range(12))
 def spec(self,terms):
  out=[set() for _ in range(12)]
  for c,word,q in terms:
   if not c:continue
   for bit,row in enumerate(self.rows(int(c),int(q))):out[bit].symmetric_difference_update(word[i] for i in self.m.bitlist(row,12))
  return [frozenset(v) for v in out]
 def product(self,left,right,targets):
  rt=[];self.m.bil_terms(rt,self.m.KT,left,right,tws_of=targets);self.E.emit_mac(rt)
 def sparse_middle(self,A,q,low,high):
  """high[i] ^= sum_j (A^q)[len(low)+i-j] low[j]."""
  m=self.m;n=len(low)+len(high)-1;pts,I=m.vinv(n)
  for a,point in enumerate(pts):
   left=self.spec([(I[j-1][a],A[j//q],q) for j in range(q,n+1,q)])
   right=self.spec([(m.gpow(point,j),word,1) for j,word in enumerate(reversed(low))])
   powers=[m.gpow(point,i) for i in range(len(high))]
   def targets(tb):return [high[i][k] for i,power in enumerate(powers) for k in m.bitlist(m.gmul(power,m.tmask(tb)),12)]
   self.product(left,right,targets)
 def triangular(self,A,q,target,inverse=False):
  if len(target)<=1:return
  cut=len(target)//2;lo,hi=target[:cut],target[cut:]
  if inverse:self.triangular(A,q,lo,True);self.sparse_middle(A,q,lo,hi);self.triangular(A,q,hi,True)
  else:self.triangular(A,q,hi);self.sparse_middle(A,q,lo,hi);self.triangular(A,q,lo)
 def self_middle(self,A,q,B):
  m=self.m;start=len(A);n=start+len(B);limit=min(len(A)-1,(n-2)//q)
  pairs=[(i,j,q*j+i-start) for j in range(1,limit+1) for i in range(1,len(A)) if start<=q*j+i<n]
  degree=(len(A)-1)+q*limit
  if len(pairs)<degree+1:
   for i,j,k in pairs:self.product(self.spec([(1,A[i],1)]),self.spec([(1,A[j],q)]),lambda tb,k=k:[B[k][b] for b in tb])
  elif pairs:
   pts,I=m.vinv(degree+1)
   for a,point in enumerate(pts):
    left=self.spec([(m.gpow(point,i),A[i],1) for i in range(1,len(A))])
    right=self.spec([(m.gpow(point,q*j),A[j],q) for j in range(1,limit+1)])
    def targets(tb):return [B[k][b] for k in range(len(B)) if start+k<=degree for b in m.bitlist(m.gmul(I[start+k][a],m.tmask(tb)),12)]
    self.product(left,right,targets)
  # The missing A[0]=1 terms are linear Frobenius maps into the high band.
  for j in range(1,len(A)):
   if tuple(A[j]) in getattr(self,'known_zero',()):continue
   index=q*j-start
   if 0<=index<len(B):
    terms=[];m.lin_rows_terms(terms,self.rows(1,q),A[j],B[index]);self.E.emit_mac(terms)
 def power(self,F,q,inverse=False):
  if len(F)<=q:return
  cut=(len(F)+1)//2;A,B=F[:cut],F[cut:]
  if inverse:
   self.power(A,q,True);self.self_middle(A,q,B);self.triangular(A,q,B,True)
  else:
   self.triangular(A,q,B);self.self_middle(A,q,B);self.power(A,q)
 def inverse(self,F):
  exponent=3;q=2
  while True:
   self.power(F,q)
   if exponent+1>=len(F):break
   q=q*q;exponent*=1+q

def build(n,optimized=False):
 m=load();b=m.Builder();F=[b.qubits(12) for _ in range(n)];sA,sB,sT=b.qubits(3);A=m.Alloc(b);E=m.Em(b,sA,sB,sT,A)
 if optimized:
  from mobius_series import BoundedArithmetic
  math=BoundedArithmetic(m,E)
 else:math=SeriesArithmetic(m,E)
 start=time.time();math.inverse(F)
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'series_inverse_{n}'+('_optimized' if optimized else '')
 tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows)
 pre=ROOT/f'parameters/{tag}.pre.kmx';kb=km.CircuitBuilder()
 for i in range(n):kb.create_quantum_register(12,f'F{i}')
 kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/'parameters/series.empty';empty.write_text('');out=ROOT/f'candidates/{tag}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
 c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb')
 result=dict(n=n,optimized=optimized,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,ops=len(rows),seconds=time.time()-start)
 print(json.dumps(result),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=35);p.add_argument('--optimized',action='store_true');a=p.parse_args();build(a.n,a.optimized)
