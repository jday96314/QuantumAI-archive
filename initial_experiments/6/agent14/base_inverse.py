"""112-Toffoli in-place GF(4096) inverse, native polynomial coordinates."""
import kickmix as km
from extension_inverse import Adapter,init_inverse12,setup
from inverse_mask_kernel import phase_inverse
from tower_inverse import map_ops,apply_map

def inplace(b,src):
 m,k,fw,bk,tr=setup();ab=Adapter(b);qq=ab.qubits(3);E=m.Em(ab,*qq,m.Alloc(ab));out=ab.qubits(12)
 init_inverse12(m,E,[q.id for q in src],out);oo=km.array([km.q(i) for i in out]);mask=b.alloc_bits(12);b.hmr(src,mask)
 apply_map(b,oo,map_ops(bk));trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));ops=map_ops(trans);apply_map(b,mask,ops)
 phase_inverse(b,oo,mask);apply_map(b,mask,ops,True);apply_map(b,oo,map_ops(bk),True)
 for x,y in zip(src,oo):b.swap(x,y)
 ab.release_recycle('base_inverse_output',out);ab.release_recycle('base_inverse_emitter',qq)
 b.bit_store(mask,False);b.free(mask)

def build():
 from fast_dqi import FastBuilder
 b=FastBuilder(km.CircuitBuilder());src=b.create_quantum_register(12,'x');inplace(b,src)
 c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
if __name__=='__main__':
 import json,numpy as np
 from pathlib import Path
 c=build();c.to_file('candidates/base_inverse12.kmb',format='kmb');m=setup()[0];x=np.arange(4096);v=np.array([m.ginv(int(a)) for a in x]);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 sim=km.Simulator(batch_size=4096,seed=142036);sim.use_same_registers_as(c);reg=c.register_data['x'];scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in set(reg)]);sim.write_across_shots(reg,bits(x));sim.do(c)
 assert np.array_equal(sim.read_across_shots(reg),bits(v));assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any();sim.do(c)
 assert np.array_equal(sim.read_across_shots(reg),bits(x));assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=4096,scope='all base-field values, zero included, inversion and inverse composition, clean scratch and exact phase');Path('results/base_inverse12_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
