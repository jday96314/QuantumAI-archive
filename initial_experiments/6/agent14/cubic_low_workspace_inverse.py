"""Existing 66-AND inverse creation with 31-AND cubic mask repayment."""
from contextlib import contextmanager
import kickmix as km
from low_workspace_inverse import norm_rt,inverse_output_rt
from cubic_mask_phase import native_phase_inverse


def bind(m):
 @contextmanager
 def inv_word(E,src):
  names=[('cubic_low_N',6),('cubic_low_M',3),('cubic_low_S',3),('cubic_low_Y',6)];wins=[E.A.open(name,n) for name,n in names];n,mm,s,y=wins;rN,rM,rS=m.norm_chain_terms(src,m.spec_rows(m.A0ROWS,src),m.spec_rows(m.A1ROWS,src),n,mm,s);rY=m.chain_rts(src,wins,src)[-2];E.macs([rN,rM,rS,rY]);E.frees([rN,rM,rS])
  for (name,_),word in reversed(list(zip(names[:3],wins[:3]))):E.A.close(name,word)
  out=E.A.open('cubic_low_out',12);E.mac(inverse_output_rt(m,src,y,out));E.free(norm_rt(m,out,y));E.A.close(names[3][0],y)
  yield out
  b=E.b.raw;mask=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in out]),mask);E.A.close('cubic_low_out',out);native_phase_inverse(b,src,mask);b.bit_store(mask,False);b.free(mask)
 return inv_word


def test():
 import json
 from pathlib import Path
 import numpy as np
 from extension_inverse import setup,Adapter
 from lean_emitter import create
 from fast_dqi import FastBuilder
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());source=b.create_quantum_register(12,'source');target=b.create_quantum_register(12,'target');ab=Adapter(b);sc,E=create(m,ab)
 with bind(m)(E,[q.id for q in source]) as iv:
  for q,r in zip(iv,target):ab.cx(q,r.id)
 ab.release_recycle('cubic_test',sc);c=b.finish_circuit();c=km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)));x=np.arange(4096);rng=np.random.default_rng(142108);y=rng.integers(0,4096,4096);want=y^np.array([m.ginv(int(v)) for v in x]);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);public=set(source)|set(target);junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public])
 for seed in [142108,142109]:
  sim=km.Simulator(batch_size=4096,seed=seed);sim.use_same_registers_as(c);sim.write_across_shots(source,bits(x));sim.write_across_shots(target,bits(y));sim.do(c);assert np.array_equal(sim.read_across_shots(source),bits(x));assert np.array_equal(sim.read_across_shots(target),bits(want));assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS',shots=8192,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like);Path('results/cubic_low_workspace_inverse_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()
