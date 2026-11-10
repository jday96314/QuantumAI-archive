"""Root-gated inverse-product parity via GF(16) cubic trace inversion.

Classical masks select a free reversible linear scaling of the divisor. The
fixed absolute-trace inverse uses a cubic Kummer tower, four-AND GF16 cubes,
and the five-product GF16 inversion core. Every quantum scratch bit is repaid.
"""
from contextlib import contextmanager
from functools import lru_cache
import json
from pathlib import Path
import numpy as np
import kickmix as km
from crt_multiplier import small_terms,indices
from tower_inverse import apply_map,map_ops,linear_xor,init_inverse,product
from fast_dqi import FastBuilder
from cubic_trace_field import data,cube_terms,inverse_network,linear


class ClassicalBuilder:
 def __init__(self,b):self.b=b
 def __getattr__(self,name):return getattr(self.b,name)
 def alloc_qubits(self,n):return self.b.alloc_bits(n)
 def reset(self,q):self.b.bit_store(q,False)
 def init_and(self,a,c,target):
  self.b.bit_store(target,False);self.b.push_condition(a);self.b.bit_store(target,True,control=c);self.b.pop_condition()
 def del_and(self,a,c,target):self.b.bit_store(target,False)


def linear_forms(b,word,mask,target):
 qs=km.array([word[j] for j in indices(mask)])
 if len(qs):b.cx(qs,target)


def xor_terms(b,left,right,target,terms):
 work=b.alloc_qubits(3);A,B,P=work
 for lm,rm,om in terms:
  linear_forms(b,left,lm,A);linear_forms(b,right,rm,B);b.init_and(A,B,P);b.cx(P,km.array([target[j] for j in indices(om)]));b.del_and(A,B,P);linear_forms(b,right,rm,B);linear_forms(b,left,lm,A)
 b.free(work)


def phase_terms(b,left,right,outcomes,terms):
 work=b.alloc_qubits(2);A,B=work
 for lm,rm,om in terms:
  linear_forms(b,left,lm,A);linear_forms(b,right,rm,B)
  for j in indices(om):b.push_condition(outcomes[j]);b.cz(A,B);b.pop_condition()
  linear_forms(b,right,rm,B);linear_forms(b,left,lm,A)
 b.free(work)


def phase_anf(b,source,outcomes,coeffs):
 for monomial,om in enumerate(coeffs):
  if not om:continue
  bits=list(indices(monomial));assert 1<=len(bits)<=2
  for j in indices(om):
   b.push_condition(outcomes[j])
   if len(bits)==1:b.z(source[bits[0]])
   else:b.cz(source[bits[0]],source[bits[1]])
   b.pop_condition()


@contextmanager
def inverse4(b,source):
 net=inverse_network();A,B=b.alloc_qubits(2);live={}
 def lin(s,target):
  if s.const:b.x(target)
  for j in s.terms:b.cx(source[j] if j<4 else live[j],target)
 for j,(a,c) in enumerate(net.nodes):
  q=b.alloc_qubits(1)[0];lin(a,A);lin(c,B);b.init_and(A,B,q);lin(c,B);lin(a,A);live[4+j]=q
 out=b.alloc_qubits(4)
 for signal,target in zip(net.outputs,out):lin(signal,target)
 yield out
 meas=b.alloc_bits(4);b.hmr(out,meas)
 for m,signal in zip(meas,net.outputs):
  b.push_condition(m)
  if signal.const:b.z(True)
  for j in signal.terms:b.z(source[j] if j<4 else live[j])
  b.pop_condition()
 b.bit_store(meas,False);b.free(meas);b.free(out)
 for j in reversed(range(len(net.nodes))):
  a,c=net.nodes[j];lin(a,A);lin(c,B);q=live.pop(4+j);b.del_and(A,B,q);lin(c,B);lin(a,A);b.free(km.array([q]))
 b.free(km.array([A,B]))


def trace_inverse(b,root,source):
 z=data();g=km.GF2Field(4);k=z['kappa'];apply_map(b,source,map_ops(z['columns']),True);a,bb,c=source[:4],source[4:8],source[8:];Q=b.alloc_qubits(4);N=b.alloc_qubits(4)
 terms=small_terms(g.modulus);qterms=[(lm,rm,g.mul(k,om)) for lm,rm,om in terms];cterms,lc,anf=cube_terms()
 xor_terms(b,bb,c,Q,qterms);linear_xor(b,a,Q,[g.mul(1<<j,1<<j) for j in range(4)])
 xor_terms(b,a,Q,N,terms)
 for word,scale in [(bb,k),(c,g.mul(k,k))]:
  xor_terms(b,word,word,N,[(lm,rm,g.mul(scale,om)) for lm,rm,om in cterms]);linear_xor(b,word,N,[g.mul(scale,v) for v in lc])
 with inverse4(b,N) as iv:
  ops=map_ops(z['trace4_columns']);apply_map(b,Q,ops);parity=b.alloc_qubits(1)[0];ladder=b.alloc_qubits(10);meas=b.alloc_bits(1)[0]
  for x,y in zip(Q,iv):b.ccx(x,y,parity)
  b.x(root);b.init_and(root[0],root[1],ladder[0])
  for j in range(1,10):b.init_and(ladder[j-1],root[j+1],ladder[j])
  b.ccz(ladder[-1],root[-1],parity)
  for j in reversed(range(1,10)):b.del_and(ladder[j-1],root[j+1],ladder[j])
  b.del_and(root[0],root[1],ladder[0]);b.x(root);b.hmr(parity,meas);b.push_condition(meas);b.cz(Q,iv);b.pop_condition();b.bit_store(meas,False);b.free(km.array([meas]));b.free(km.array([parity]));b.free(ladder);apply_map(b,Q,ops,True)
 nm=b.alloc_bits(4);b.hmr(N,nm);phase_terms(b,a,Q,nm,terms)
 for word,scale in [(bb,k),(c,g.mul(k,k))]:phase_anf(b,word,nm,[g.mul(scale,v) for v in anf])
 b.bit_store(nm,False);b.free(nm);b.free(N)
 qm=b.alloc_bits(4);b.hmr(Q,qm);phase_terms(b,bb,c,qm,qterms);phase_anf(b,a,qm,[g.mul(1<<j.bit_length()-1,1<<j.bit_length()-1) if j and j&(j-1)==0 else 0 for j in range(16)])
 b.bit_store(qm,False);b.free(qm);b.free(Q);apply_map(b,source,map_ops(z['columns']))


@lru_cache(None)
def product_columns():
 f=km.GF2Field(12)
 return tuple(tuple(f.mul(1<<i,1<<j) for i in range(12)) for j in range(12))


def classical_product_xor(b,source,coefficient,target):
 for control,cols in zip(coefficient,product_columns()):
  b.push_condition(control);linear_xor(b,source,target,cols);b.pop_condition()


def apply(b,root,a,c,mask):
 active=b.alloc_bits(1)[0];b.bit_store(active,True,control=mask);b.push_condition(active)
 ops=map_ops(data()['dual_columns']);apply_map(b,mask,ops);ik=b.alloc_bits(12);sc=init_inverse(ClassicalBuilder(b),mask,ik)
 value=b.alloc_qubits(12);product(b,a,c,value);scaled=b.alloc_qubits(12);classical_product_xor(b,value,ik,scaled);classical_product_xor(b,scaled,mask,value);b.free(value)
 trace_inverse(b,root,scaled)
 value=b.alloc_qubits(12);classical_product_xor(b,scaled,mask,value);classical_product_xor(b,value,ik,scaled);b.free(scaled);b.del_gf2_mul(a,c,target=value,field=km.GF2Field(12));b.free(value)
 b.bit_store(ik,False);b.bit_store(sc,False);b.free(ik);b.free(sc);apply_map(b,mask,ops,True);b.pop_condition();b.bit_store(active,False);b.free(km.array([active]))


def test():
 b=FastBuilder(km.CircuitBuilder());root=b.create_quantum_register(12,'root');a=b.create_quantum_register(12,'a');c=b.create_quantum_register(12,'c');mask=b.create_classical_register(12,'mask');apply(b,root,a,c,mask);circ=b.finish_circuit();n=4096;rng=np.random.default_rng(71319);av=np.arange(n);cv=rng.integers(0,n,n);rv=rng.integers(0,n,n);rv[::2]=0;mv=rng.integers(0,n,n);mv[::31]=0;bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 sim=km.Simulator(batch_size=n,seed=71);sim.use_same_registers_as(circ)
 for reg,v in [(root,rv),(a,av),(c,cv),(mask,mv)]:sim.write_across_shots(reg,bits(v))
 sim.do(circ);f=km.GF2Field(12);want=np.array([r==0 and (f.pow(f.mul(int(x),int(y)),4094)&int(m)).bit_count()%2 for r,x,y,m in zip(rv,av,cv,mv)],dtype=bool);bad=np.flatnonzero(sim.read_phase_flipped_across_shots()!=want);failures={}
 for name,reg,v in [('root',root,rv),('a',a,av),('c',c,cv),('mask',mask,mv)]:
  errors=np.flatnonzero((sim.read_across_shots(reg)!=bits(v)).any(0))
  if len(errors):failures[name]=errors[:8].tolist()
 public=set(root)|set(a)|set(c);scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));out=dict(status='PASS' if not len(bad) and not failures and not len(junk) else 'FAIL',scope='all first-operand field values, random second operand and runtime masks, zero/nonzero roots, zero masks, exact phase, all input/scratch restoration',shots=n,Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like,phase_errors=bad[:8].tolist(),failures=failures,scratch_errors=junk[:8].tolist());Path('results/cubic_inverse_product_kernel_check.json').write_text(json.dumps(out,indent=2)+'\n');circ.to_file('candidates/cubic_inverse_product_kernel.kmb',format='kmb');print(json.dumps(out),flush=True)
if __name__=='__main__':test()
