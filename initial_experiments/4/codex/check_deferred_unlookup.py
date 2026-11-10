"""Check joint ANF repair with lookup outcomes measured and reused earlier."""
import json
import numpy as np
import kickmix as km
from monomial_unlookup import erase

rng=np.random.default_rng(717219);n=9;width=46;extra_width=4;N=1024
table=rng.integers(0,1<<(width+extra_width),1<<n,dtype=np.int64).tolist()
b=km.CircuitBuilder();a=b.create_quantum_register(n,'a');target=b.create_quantum_register(width,'target');early=b.alloc_qubits(extra_width)
outcomes=b.alloc_bits(extra_width);b.hmr(early,outcomes)
# The measured qubits are genuinely reusable before final phase repair.
b.init_and(a[0],a[1],early[0]);b.del_and(a[0],a[1],early[0])
erase(b,a,target,table,extra_outcomes=outcomes);b.free(early);c=b.finish_circuit()
addresses=np.arange(N)%(1<<n);values=np.array(table,dtype=np.int64)[addresses]
bits=lambda x,w:((x[None,:]>>np.arange(w)[:,None])&1).astype(bool)
s=km.Simulator(batch_size=N,seed=3361);s.use_same_registers_as(c);s.write_across_shots(a,bits(addresses,n));s.write_across_shots(target,bits(values,width));s.write_across_shots(early,bits(values>>width,extra_width));s.do(c)
assert np.array_equal(s.read_across_shots(a),bits(addresses,n));assert not s.read_across_shots(target).any();assert not s.read_phase_flipped_across_shots().any()
public=set(a)|set(target);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not s.read_across_shots(scratch).any()
result=dict(status='PASS',shots=N,addresses=1<<n,late_bits=width,early_bits=extra_width,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)
open('results/deferred_unlookup_check.json','w').write(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
