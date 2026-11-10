"""Add a short signed carry whose sign is known from the coordinate."""
import kickmix as km

def add(b,counter,target,negative,subtract=False,qubit_limit=None,workspace=()):
 assert 1<=len(counter)<=2
 # The modular carry is nonpositive exactly on the known negative branch.
 # Negating a two-bit number preserves bit 0 and toggles bit 1 on bit 0.
 if len(counter)==2:b.ccx(negative,counter[0],counter[1])
 b.cx(negative,target)
 from short_adder import add as unsigned_add
 unsigned_add(b,counter,target,subtract=subtract,qubit_limit=qubit_limit,carry_workspace=workspace)
 b.cx(negative,target)
 if len(counter)==2:b.ccx(negative,counter[0],counter[1])

def test():
 import json,numpy as np
 rng=np.random.default_rng(2274);results=[]
 for width in [1,2]:
  for sub in [False,True]:
   b=km.CircuitBuilder();s=b.create_quantum_register(1,'s')[0];c=b.create_quantum_register(width,'c');t=b.create_quantum_register(10,'t');add(b,c,t,s,sub);circ=b.finish_circuit();N=2048
   sign=rng.integers(0,2,N);mag=rng.integers(0,1<<width,N);delta=np.where(sign,-mag,mag);cc=delta%(1<<width);tt=rng.integers(0,1024,N)
   bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
   sim=km.Simulator(batch_size=N,seed=915);sim.use_same_registers_as(circ)
   for reg,v in [(km.array([s]),sign),(c,cc),(t,tt)]:sim.write_across_shots(reg,bits(v,len(reg)))
   sim.do(circ);expected=(tt-delta if sub else tt+delta)%1024
   for reg,v in [(km.array([s]),sign),(c,cc),(t,expected)]:assert np.array_equal(sim.read_across_shots(reg),bits(v,len(reg)))
   assert not sim.read_phase_flipped_across_shots().any();public=set(c)|set(t)|{s};scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(scratch).any()
   results.append(dict(width=width,subtract=sub,shots=N,Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like))
 out=dict(status='PASS',cases=results);open('results/signed_counter_check.json','w').write(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':test()
