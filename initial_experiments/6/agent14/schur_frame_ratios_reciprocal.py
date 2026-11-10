"""Recreate the extraction frame directly from affine head ratios."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from tower_inverse import apply_map,map_ops
from fast_dqi import FastBuilder
from spooky_extension_norm_inverse import init_delayed,inverse_mask_native


def build(N):
 m,k,fw,bk,tr=setup();m.inv_word=bind(m);n=N//2;b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(N)];g=b.create_quantum_register(12,'g');ghost=b.create_classical_register(12,'frame_ghost');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));G=[q.id for q in g];even=b.alloc_qubits(1)[0]
 def parity():
  if not n%2:b.x(even)
  for i in range(1,n):b.flip_if_equal(km.array(list(S[2*i])+list(S[2*i+1])),1<<12,target=even)
 def swaps(back=False):
  for i in (reversed(range(0,n-1,2)) if back else range(0,n-1,2)):
   for q,r in zip(S[2*i+1],S[2*i+3]):b.cswap(even,q,r)
 parity();swaps();b.x(g[0])
 for i in range(0,n,2):
  word=[q.id for q in S[2*i+1]]
  with m.inv_word(E,word) as wi,m.opened(E,'frame_ratio_product',12) as tmp:m.mon_slot(E,G,word,wi,tmp)
 swaps(True)
 out=init_delayed(m,E,G)
 b.x(even)
 for q,r in zip(g,out):b.cswap(even,q,km.q(r))
 b.x(even)
 sq=map_ops([m.gmul(1<<i,1<<i) for i in range(12)]);outreg=km.array([km.q(q) for q in out]);apply_map(b,outreg,sq)
 for q,bit in zip(out,ghost):ab.z(q,bit)
 b.bit_store(ghost,False);apply_map(b,outreg,sq,True)
 mask=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in out]),mask);E.A.close('ext_delayed_inverse',out);inverse_mask_native(b,G,mask,fw,bk);b.bit_store(mask,False);b.free(mask)
 parity();b.free(km.array([even]));apply_map(b,g,map_ops([m.gmul(1<<i,1<<i) for i in range(12)]));ab.release_recycle('ratio_frame_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(N):
 m=setup()[0];shots=4096;n=N//2;rng=np.random.default_rng(142095+N);values=rng.integers(0,4096,(shots,N));want=[]
 for row,v in enumerate(values):
  d=row%(n+1);raw=rng.integers(1,4096,d).tolist();g=1
  for i in range(n):g=m.gmul(m.gmul(raw[i],raw[i]) if i<d else 1,m.ginv(g))
  v[1]=raw[0] if d else 1
  for i in range(1,n):
   if i<d:v[2*i+1]=m.gmul(raw[i],m.ginv(raw[i-1]))
   else:v[2*i]=0;v[2*i+1]=1
  want.append(m.ginv(g))
 c=build(N);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142095);sim.use_same_registers_as(c)
 for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(values[:,i]))
 sim.do(c);assert np.array_equal(sim.read_across_shots(c.register_data['g']),bits(want))
 for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(values[:,i]))
 public={q for name,reg in c.register_data.items() if name!='frame_ghost' for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any();r=dict(status='PASS',N=N,shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,scope='all completed-prefix lengths, arbitrary active affine heads and ratios, restored source, reciprocal frame with regular-frame ghost payment, zero input ghost in standalone test');Path(f'results/schur_frame_ratios_reciprocal_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 for N in [6,70]:test(N)
