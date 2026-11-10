"""Unequal-width controlled additions using dirty input bits as zero padding.
Subtract their restored high contribution after the main Cuccaro sweep.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from independent_arithmetic import controlled_add

def main():
 rows=[];rng=np.random.default_rng(725)
 for m,n in [(1,7),(3,7),(6,7),(5,19),(17,19)]:
  for carry in [False,True]:
   for sub in [False,True]:
    b=km.CircuitBuilder();a=b.create_quantum_register(m,'a');t=b.create_quantum_register(n,'t');d=b.create_quantum_register(n-m,'d');ctrl=b.create_quantum_register(1,'ctrl')[0];ci=b.create_quantum_register(1,'ci')[0] if carry else None
    controlled_add(b,a,t,ctrl,sub,ci,d);c=b.finish_circuit();N=4096
    av=rng.integers(0,1<<m,N);tv=rng.integers(0,1<<n,N);dv=rng.integers(0,1<<(n-m),N);cv=rng.integers(0,2,N);iv=rng.integers(0,2,N) if carry else np.zeros(N,dtype=int)
    bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool)
    s=km.Simulator(batch_size=N,seed=12);s.use_same_registers_as(c)
    for reg,v in [(a,av),(t,tv),(d,dv),(km.array([ctrl]),cv)]+([(km.array([ci]),iv)] if carry else []):s.write_across_shots(reg,bits(v,len(reg)))
    s.do(c);expected=(tv+(-1 if sub else 1)*cv*(av+iv))%(1<<n)
    for reg,v in [(a,av),(t,expected),(d,dv),(km.array([ctrl]),cv)]+([(km.array([ci]),iv)] if carry else []):assert np.array_equal(s.read_across_shots(reg),bits(v,len(reg)))
    assert not s.read_phase_flipped_across_shots().any()
    public=set(a)|set(t)|set(d)|{ctrl}|({ci} if carry else set());scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not s.read_across_shots(scratch).any()
    rows.append(dict(m=m,n=n,carry=carry,subtract=sub,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=N))
 Path('results/dirty_borrow_tests.json').write_text(json.dumps(dict(status='PASS',rows=rows),indent=2)+'\n');print(json.dumps(dict(status='PASS',cases=len(rows),shots=len(rows)*4096)))
if __name__=='__main__':main()
