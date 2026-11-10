"""Exhaustive field/basis-mask and exhaustive classical-mask checks."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from cubic_inverse_product_kernel import apply
from fast_dqi import FastBuilder

b=FastBuilder(km.CircuitBuilder());root=b.create_quantum_register(12,'root');a=b.create_quantum_register(12,'a');c=b.create_quantum_register(12,'c');mask=b.create_classical_register(12,'mask');apply(b,root,a,c,mask);circ=b.finish_circuit();rng=np.random.default_rng(711937);f=km.GF2Field(12);values=np.arange(4096);public=set(root)|set(a)|set(c);scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public]);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
cases=[(f'basis_mask_{j}',values,np.ones(4096,dtype=int),np.zeros(4096,dtype=int),np.full(4096,1<<j)) for j in range(12)]
cases += [('all_masks',rng.integers(0,4096,4096),rng.integers(0,4096,4096),np.zeros(4096,dtype=int),values),('all_zero_divisors',np.zeros(4096,dtype=int),values,np.zeros(4096,dtype=int),values),('nonroots',values,rng.integers(0,4096,4096),rng.integers(1,4096,4096),values)]
for seed,(name,av,cv,rv,mv) in enumerate(cases):
 sim=km.Simulator(batch_size=4096,seed=719+seed);sim.use_same_registers_as(circ)
 for reg,v in [(root,rv),(a,av),(c,cv),(mask,mv)]:sim.write_across_shots(reg,bits(v))
 sim.do(circ);want=np.array([r==0 and (f.pow(f.mul(int(x),int(y)),4094)&int(m)).bit_count()%2 for r,x,y,m in zip(rv,av,cv,mv)],dtype=bool)
 assert np.array_equal(sim.read_phase_flipped_across_shots(),want),name
 for reg,v in [(root,rv),(a,av),(c,cv),(mask,mv)]:assert np.array_equal(sim.read_across_shots(reg),bits(v)),name
 assert not sim.read_across_shots(scratch).any(),name
 reports.append(dict(name=name,status='PASS',shots=4096))
result=dict(status='PASS',scope='all4096 divisors against all12 mask basis vectors; all4096 runtime masks; zero divisors and nonroot paths; source/classical-mask/scratch restoration and exact phase',shots=4096*len(cases),Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like,cases=reports);Path('results/cubic_inverse_product_kernel_exhaustive.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
