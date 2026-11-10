"""Original low-width controlled Cuccaro arithmetic and native kernel checks.
No incumbent source imported. MAJ followed by conditional UMA restores source.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km

def controlled_add(b,a,t,ctrl,subtract=False,carry_in=None,borrow=()):
 """t += ctrl*(a+carry_in), preserving a and carry_in; m <= n.

Padding implements unequal widths. Counts include all executed Toffolis.
"""
 assert 0<len(a)<=len(t)
 if subtract:b.x(t)
 need=len(t)-len(a);borrowed=km.array(list(borrow)[:need]);pad=b.alloc_qubits(need-len(borrowed));padding=km.array([*borrowed,*pad]);aa=km.array([*a,*padding])
 own=b.alloc_qubits(int(carry_in is None));initial=own[0] if len(own) else carry_in
 for i in range(len(t)-1):
  c=initial if i==0 else aa[i-1]
  b.cx(aa[i],t[i]);b.cx(aa[i],c);b.ccx(c,t[i],aa[i])
 # Conditional high sum, shared control reduces two ANDs to one.
 top=aa[-1];boundary=initial if len(t)==1 else aa[-2]
 b.cx(boundary,top);b.ccx(ctrl,top,t[-1]);b.cx(boundary,top)
 for i in reversed(range(len(t)-1)):
  c=initial if i==0 else aa[i-1]
  b.ccx(c,t[i],aa[i]);b.cx(aa[i],c);b.cx(aa[i],t[i])
  b.cx(aa[i],c);b.ccx(ctrl,c,t[i]);b.cx(aa[i],c)
 b.free(own)
 if len(borrowed):controlled_add(b,borrowed,t[len(a):],ctrl,True)
 b.free(pad)
 if subtract:b.x(t)

def rounded_product(b,a,z,t,exponent,subtract=False,borrow=()):
 """Accumulate signed two's complement z times a, with each row rounded."""
 for j,ctrl in enumerate(z):
  e=exponent+j;sign=subtract^(j==len(z)-1)
  lo=max(0,e);src=a[max(0,-e):min(len(a),len(t)-e)]
  rounding=a[-e-1] if e<0 and -e-1<len(a) else None
  available=[q for q in borrow if q!=ctrl and q not in set(a) and q not in set(t)]
  if len(src):controlled_add(b,src,t[lo:],ctrl,sign,rounding,available)
  elif rounding is not None:controlled_add(b,km.array([rounding]),t,ctrl,sign,borrow=available)

def tests():
 rows=[];rng=np.random.default_rng(9027)
 for m,n in [(1,1),(1,7),(3,7),(7,7),(17,20)]:
  for carry in [False,True]:
   for sub in [False,True]:
    b=km.CircuitBuilder();a=b.create_quantum_register(m,'a');t=b.create_quantum_register(n,'t');ctrl=b.create_quantum_register(1,'ctrl')[0];ci=b.create_quantum_register(1,'ci')[0] if carry else None
    controlled_add(b,a,t,ctrl,sub,ci);c=b.finish_circuit();N=4096
    av=rng.integers(0,1<<m,N);tv=rng.integers(0,1<<n,N);cv=rng.integers(0,2,N);iv=rng.integers(0,2,N) if carry else np.zeros(N,dtype=int)
    if m==n==1:
     a0,t0,c0,i0=np.indices((2,2,2,2)).reshape(4,-1);av=np.resize(a0,N);tv=np.resize(t0,N);cv=np.resize(c0,N);iv=np.resize(i0,N) if carry else np.zeros(N,dtype=int)
    bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool)
    s=km.Simulator(batch_size=N,seed=17,counting_operations=True);s.use_same_registers_as(c)
    s.write_across_shots(a,bits(av,m));s.write_across_shots(t,bits(tv,n));s.write_across_shots(km.array([ctrl]),bits(cv,1))
    if carry:s.write_across_shots(km.array([ci]),bits(iv,1))
    s.do(c);expected=(tv+(-1 if sub else 1)*cv*(av+iv))% (1<<n)
    assert np.array_equal(s.read_across_shots(t),bits(expected,n)),(m,n,carry,sub)
    assert np.array_equal(s.read_across_shots(a),bits(av,m))
    assert not s.read_phase_flipped_across_shots().any()
    public=set(a)|set(t)|{ctrl}|({ci} if carry else set());scratch=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);assert not s.read_across_shots(scratch).any()
    rows.append(dict(m=m,n=n,carry=carry,subtract=sub,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=N))
 Path('results/controlled_cuccaro_tests.json').write_text(json.dumps(dict(status='PASS',rows=rows),indent=2)+'\n');print(json.dumps(dict(status='PASS',cases=len(rows),native_shots=len(rows)*4096)))
if __name__=='__main__':tests()
