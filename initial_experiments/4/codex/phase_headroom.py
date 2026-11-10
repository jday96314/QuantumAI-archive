"""Phase a shifted coefficient using conditionally clean normalized-input bits.

The first h input bits are zero whenever the normalization shift is at least h.
Barrel stage 2^j touches at most 2^(j+1)-1 headroom bits, which are zero whenever
that stage executes. Headroom above the actual shift can remain dirty; threshold
controls exclude it from the phase.
"""
from fractions import Fraction
import kickmix as km

def apply(b,y,core,shift):
 f=len(y)
 for j in range(core):b.z_pow(y[j],Fraction(2)**(j-f+1))
 for h,q in enumerate(y[core:]):
  threshold=h+1;low=(threshold&-threshold).bit_length()-1
  assert threshold<(1<<len(shift))
  current=frozenset([shift[low]]);A=b.alloc_qubits(1)[0];records=[]
  def parity(nodes):
   for node in nodes:b.cx(node,A)
  for j in range(low+1,len(shift)):
   t=b.alloc_qubits(1)[0];parity(current);b.init_and(shift[j],A,t);parity(current)
   records.append((current,shift[j],t))
   current=frozenset([t]) if threshold>>j&1 else current^frozenset([shift[j],t])
  parity(current);angle=Fraction(2)**(core+h-f+1)
  if angle==1:b.cz(A,q)
  else:
   t=b.alloc_qubits(1)[0];b.init_and(A,q,t);b.z_pow(t,angle);b.del_and(A,q,t);b.free(km.array([t]))
  parity(current)
  for old,s,t in reversed(records):parity(old);b.del_and(s,A,t);parity(old);b.free(km.array([t]))
  b.free(km.array([A]))

def test():
 import json,numpy as np
 from binary_centered_invsqrt import cshift
 b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');y0=b.create_quantum_register(26,'y');flags=b.create_quantum_register(5,'flags')
 y=km.array([*y0,x[21],*x[:10]]);occ=[27,28,30,34]
 for k in range(4):cshift(b,y,1<<k,flags[k+1],occ[k])
 apply(b,y,27,flags[1:])
 for k in reversed(range(4)):cshift(b,y,1<<k,flags[k+1],occ[k],True)
 c=b.finish_circuit();N=2048;rng=np.random.default_rng(71881);kk=np.arange(N)%22
 xx=(rng.integers(0,1<<21,N)|(1<<21))&~((1<<kk)-1);yy=rng.integers(0,1<<26,N)
 bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
 sim=km.Simulator(batch_size=N,seed=819);sim.use_same_registers_as(c)
 for reg,values in [(x,xx),(y0,yy),(flags,kk)]:sim.write_across_shots(reg,bits(values,len(reg)))
 sim.do(c)
 for reg,values in [(x,xx),(y0,yy),(flags,kk)]:assert np.array_equal(sim.read_across_shots(reg),bits(values,len(reg)))
 expected=((yy+(1<<26))<<(kk//2)).astype(np.uint64)<<(64-37)
 actual=np.array([int(sim.read_shot_phase(i)*(1<<63))&((1<<64)-1) for i in range(N)],dtype=np.uint64)
 assert np.array_equal(actual,expected),(actual[:8],expected[:8])
 public=set(x)|set(y0)|set(flags);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(scratch).any()
 out=dict(status='PASS',shots=N,normalization_shifts=22,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)
 print(json.dumps(out),flush=True);open('results/phase_headroom_check.json','w').write(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':test()
