"""Trade controlled-increment gates for fewer clean carry wires.

A low block's carry is computed, its intermediate ANDs measured away, and the
high block incremented. Recomputing the low prefix then erases that carry.
Only then does the low block change. Every erasure has its ordinary CZ repair.
"""
from functools import lru_cache
import kickmix as km
import short_adder
ORIGINAL=short_adder.add

@lru_cache(None)
def plan(n,available):
 if n<=2 or available>=n-2:return (max(0,n-1),0)
 if available<=0:return (float('inf'),0)
 best=(float('inf'),0)
 for k in range(1,min(n-1,available)+1):
  cost=2*k-1+plan(n-k,available-1)[0]+plan(k,available)[0]
  if cost<best[0]:best=(cost,k)
 return best

def increment(b,c,target,workspace,limit):
 n=len(target);available=len(workspace)+limit-b.num_allocated_qubits
 cost,k=plan(n,available)
 if cost==float('inf'):raise ValueError(f'increment needs more workspace: n={n}, available={available}')
 if not k:return ORIGINAL(b,km.array([c]),target,carry_workspace=workspace)
 take=min(k,len(workspace));own=b.alloc_qubits(k-take);chain=km.array([*workspace[:take],*own]);flag=chain[-1]
 for i in range(k):b.init_and(c if i==0 else chain[i-1],target[i],chain[i])
 for i in reversed(range(k-1)):b.del_and(c if i==0 else chain[i-1],target[i],chain[i])
 # Release owned prefix wires, retaining only the carry itself.
 owned_flag=flag in set(own);discard=km.array([q for q in own if q!=flag]);b.free(discard)
 remaining=km.array([q for q in workspace if q!=flag])
 increment(b,flag,target[k:],remaining,limit)
 # Low target bits have not changed, so the carry can now be erased exactly.
 count=k-1;take2=min(count,len(remaining));own2=b.alloc_qubits(count-take2);prefix=km.array([*remaining[:take2],*own2])
 for i in range(count):b.init_and(c if i==0 else prefix[i-1],target[i],prefix[i])
 b.del_and(c if not count else prefix[-1],target[k-1],flag)
 for i in reversed(range(count)):b.del_and(c if i==0 else prefix[i-1],target[i],prefix[i])
 b.free(own2)
 if owned_flag:b.free(km.array([flag]))
 increment(b,c,target[:k],workspace,limit)

def add(b,lhs,target,subtract=False,carry_in=None,carry_workspace=(),prefix_bits=0,qubit_limit=None):
 if qubit_limit is not None and len(lhs)==1 and carry_in is None:
  available=len(carry_workspace)+qubit_limit-b.num_allocated_qubits
  if available<len(target)-2:
   if plan(len(target),available)[0]==float('inf'):raise ValueError('increment cannot fit requested qubit limit')
   if subtract:b.x(target)
   increment(b,lhs[0],target,km.array(carry_workspace),qubit_limit)
   if subtract:b.x(target)
   return
 if qubit_limit is not None and len(lhs)>1:
  m,n=len(lhs),len(target);fits=False
  for candidate in range(min(m,n-1)+1):
   if candidate==0:extra=max(0,n-2-len(carry_workspace))
   else:
    initial=int(carry_in is None and not len(carry_workspace));remaining=len(carry_workspace)-int(carry_in is None and bool(len(carry_workspace)))
    extra=initial+max(0,n-candidate-2-remaining)
   fits|=b.num_allocated_qubits+extra<=qubit_limit
  if not fits:
   k=min(m,n-1);own=km.array([]);available=km.array(carry_workspace)
   if carry_in is not None:initial=carry_in
   elif len(available):initial=available[0];available=available[1:]
   else:own=b.alloc_qubits(1);initial=own[0]
   room=len(available)+qubit_limit-b.num_allocated_qubits
   if room<0:
    if len(own):b.free(own)
    raise ValueError('prefix initial carry cannot fit requested qubit limit')
   if m==k and plan(n-k,room)[0]==float('inf'):
    if len(own):b.free(own)
    raise ValueError('full prefix cannot fit requested qubit limit')
   if subtract:b.x(target)
   for i in range(k):
    c=initial if i==0 else lhs[i-1];b.cx(lhs[i],target[i]);b.cx(lhs[i],c);b.ccx(c,target[i],lhs[i])
   boundary=lhs[k-1]
   if m>k:b.cx(lhs[k],target[k]);b.cx(boundary,target[k])
   else:increment(b,boundary,target[k:],available,qubit_limit)
   for i in reversed(range(k)):
    c=initial if i==0 else lhs[i-1];b.ccx(c,target[i],lhs[i]);b.cx(lhs[i],c);b.cx(c,target[i])
   if len(own):b.free(own)
   if subtract:b.x(target)
   return
 return ORIGINAL(b,lhs,target,subtract,carry_in,carry_workspace,prefix_bits,qubit_limit)

def install():
 import booth_adder
 short_adder.add=add;booth_adder.add=add

def test():
 import json,numpy as np
 rng=np.random.default_rng(6119);results=[]
 for n,d in [(5,2),(7,3),(10,4),(12,4),(18,6),(24,7)]:
  for subtract in [False,True]:
   b=km.CircuitBuilder();c=b.create_quantum_register(1,'c');t=b.create_quantum_register(n,'t');limit=n+1+d
   add(b,c,t,subtract=subtract,qubit_limit=limit);circuit=b.finish_circuit();N=2048
   cc=rng.integers(0,2,N);tt=rng.integers(0,1<<n,N);tt[:4]=[0,1,(1<<n)-1,(1<<n)-2]
   bits=lambda x,w:((x[None,:]>>np.arange(w)[:,None])&1).astype(bool)
   sim=km.Simulator(batch_size=N,seed=41);sim.use_same_registers_as(circuit);sim.write_across_shots(c,bits(cc,1));sim.write_across_shots(t,bits(tt,n));sim.do(circuit)
   expected=(tt-cc if subtract else tt+cc)%(1<<n)
   assert np.array_equal(sim.read_across_shots(t),bits(expected,n));assert np.array_equal(sim.read_across_shots(c),bits(cc,1));assert not sim.read_phase_flipped_across_shots().any()
   scratch=km.array([km.q(i) for i in range(n+1,circuit.num_qubits)]);assert not sim.read_across_shots(scratch).any();assert circuit.num_qubits<=limit
   results.append(dict(n=n,available=d,subtract=subtract,Q=circuit.num_qubits,T=circuit.operation_statistics().num_toffoli_like,shots=N))
 out=dict(status='PASS',cases=results);print(json.dumps(out),flush=True);open('results/space_increment_check.json','w').write(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':test()
