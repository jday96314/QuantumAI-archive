"""Borrow normalized input bits for the rounded square coordinate.

The normalized top bit starts and ends at one. Its temporary zero state is
used for the rounded endpoint guard. HMR erasure streams a 9-bit ANF bank.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from short_adder import add
from rounded_square import apply as square
from streamed_unlookup import erase


def apply(b,x,h,inverse=False):
 assert len(x)==22 and len(h)==11
 def fold():b.x(x[:13]);b.cx(x[13],x[:13])
 fold();b.x(x[21]);u=km.array(list(x[5:13])+[x[21]]);add(b,x[4:5],u)
 if inverse:
  extra=b.alloc_qubits(4);target=km.array(list(h)+list(extra));table=[((i*i+32)//64)&2047 for i in range(512)];erase(b,u,target,table,split=4);b.free(extra)
 else:square(b,u,h)
 add(b,x[4:5],u,subtract=True);b.x(x[21]);fold()


def build(inverse=False):
 b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');h=b.create_quantum_register(11,'h');apply(b,x,h,inverse);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test():
 rng=np.random.default_rng(142080);lo=np.arange(16384);xx=lo|(rng.integers(0,128,len(lo))<<14)|(1<<21);z=lo-8192;u=np.abs((z+16)>>5);want=(u*u+32)>>6;bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);cf=build();ci=build(True);rows=[]
 for seed in (142080,142081):
  for inverse,c in [(False,cf),(True,ci)]:
   sim=km.Simulator(batch_size=len(xx),seed=seed);sim.use_same_registers_as(c);regs=c.register_data;sim.write_across_shots(regs['x'],bits(xx,22))
   if inverse:sim.write_across_shots(regs['h'],bits(want,11))
   sim.do(c);bad=np.flatnonzero((sim.read_across_shots(regs['h'])!=bits(np.zeros_like(want) if inverse else want,11)).any(0));assert np.array_equal(sim.read_across_shots(regs['x']),bits(xx,22));public=set(regs['x'])|set(regs['h']);sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(seed=seed,inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=bad[:8].tolist(),junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not len(bad) and not len(junk) and not len(phase)
   sim.do(cf if inverse else ci);assert np.array_equal(sim.read_across_shots(regs['x']),bits(xx,22));assert np.array_equal(sim.read_across_shots(regs['h']),bits(want if inverse else np.zeros_like(want),11));assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',shots=4*len(xx),rows=rows,scope='all 16384 centered coordinates, random high mantissa bits, normalized top=one, rounded endpoint256, both directions and compositions, input/scratch/phase restoration');Path('results/rounded_coordinate_inplace_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':test()
