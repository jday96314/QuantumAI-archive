"""Exact rounded absolute centered coordinate and narrow square from mantissa."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from short_adder import add
from rounded_square import apply as square


def coordinate(b,x,u,inverse=False):
 def fold():
  # sign(z)=not x13, z=(M mod 2^14)-8192. One's-complement
  # folding yields -z-1 on the negative branch and z otherwise.
  b.x(x[:13]);b.cx(x[13],x[:13])
 fold()
 if inverse:add(b,x[4:5],u,subtract=True);b.cx(x[5:13],u[:8])
 else:b.cx(x[5:13],u[:8]);add(b,x[4:5],u)
 fold()


def apply(b,x,h,inverse=False):
 u=b.alloc_qubits(9);coordinate(b,x,u);square(b,u,h,inverse);coordinate(b,x,u,True);b.free(u)


def build(inverse=False):
 b=km.CircuitBuilder();x=b.create_quantum_register(14,'x');h=b.create_quantum_register(11,'h');apply(b,x,h,inverse);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test():
 cf=build();ci=build(True);xx=np.arange(16384);z=xx-8192;u=np.abs((z+16)>>5);want=(u*u+32)>>6;bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);rows=[]
 for rev,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=len(xx),seed=142056+rev);sim.use_same_registers_as(c);regs=c.register_data;sim.write_across_shots(regs['x'],bits(xx,14))
  if rev:sim.write_across_shots(regs['h'],bits(want,11))
  sim.do(c);bad=np.flatnonzero((sim.read_across_shots(regs['h'])!=bits(np.zeros_like(want) if rev else want,11)).any(0));assert np.array_equal(sim.read_across_shots(regs['x']),bits(xx,14));public=set(regs['x'])|set(regs['h']);sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(inverse=rev,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=bad[:8].tolist(),junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not len(bad) and not len(junk) and not len(phase)
  sim.do(cf if rev else ci);assert np.array_equal(sim.read_across_shots(regs['x']),bits(xx,14));assert np.array_equal(sim.read_across_shots(regs['h']),bits(want if rev else np.zeros_like(want),11));assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',shots=32768,rows=rows,scope='all 16384 centered coordinates, both endpoints, rounded absolute-square value, fresh inverse and both compositions, restored input and clean scratch/phase');Path('results/rounded_coordinate_native.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':test()
