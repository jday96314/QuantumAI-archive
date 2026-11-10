"""Pair an inverse with delayed repair of its old-input measurements.

Forward HMR outcomes stay private and live until the matching reverse
inversion reconstructs their controls. The intervening computation must
be a basis permutation/diagonal phase. No basis value is measured in Z.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup,init_inverse12,inverse_phase,extension_norm_rt


def build(reverse=False):
 m,k,fw,bk,tr=setup();b=FastBuilder(km.CircuitBuilder());x=b.alloc_qubits(24);headers=''.join(f'APPEND_TO_REGISTER q{q.id} r1\n' for q in x)+'REGISTER r1 "x"\n';ghost=b.create_classical_register(24,'ghost');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));src=[q.id for q in x]
 with m.opened(E,'extension_norm',12) as nn,m.opened(E,'extension_norm_inverse',12) as ni:
  nr=extension_norm_rt(m,src,nn,k);E.mac(nr);init_inverse12(m,E,nn,ni);out=ab.qubits(24)
  left=[a^c for a,c in zip(m.spec_plain(src[:12]),m.spec_plain(src[12:]))];rt=[];m.bil_terms(rt,m.KT,left,m.spec_plain(ni),out[:12]);E.mac(rt);E.mac(m.prod_rt(out[12:],src[12:],ni));E.free(extension_norm_rt(m,out,ni,k));E.free(nr)
 if reverse:
  for q,mask in zip(out,ghost):ab.z(q,mask)
  measurement=b.alloc_bits(24);b.hmr(x,measurement);b.free(x);inverse_phase(b,E,out,measurement,k);b.bit_store(measurement,False);b.free(measurement);b.bit_store(ghost,False);reclaimed=b.alloc_qubits(24);assert {q.id for q in reclaimed}==set(src),([q.id for q in reclaimed],src)
 else:b.bit_store(ghost,False);b.hmr(x,ghost)
 for q,r in zip(x,out):ab.swap(q.id,r)
 ab.release_recycle('spooky_output',out);ab.release_recycle('spooky_emitter',sc);c=b.finish_circuit();held={v.id for v in ghost}
 return km.Circuit(headers+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits) if i not in held))

def test():
 cf=build();ci=build(True);m,k,fw,bk,tr=setup();rng=np.random.default_rng(142046);shots=16384;vals=rng.integers(0,1<<24,shots);vals[:4096]=np.arange(4096);vals[4096:8192]=np.arange(4096)<<12
 def inv(v):
  a,b=int(v)&4095,int(v)>>12;n=m.gmul(a,a)^m.gmul(a,b)^m.gmul(k,m.gmul(b,b));ni=m.ginv(n)
  return m.gmul(a^b,ni)|(m.gmul(b,ni)<<12)
 want=np.array([inv(v) for v in vals]);bits=lambda v:((v[None,:]>>np.arange(24)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142046);sim.use_same_registers_as(ci);x=cf.register_data['x'];ghost=cf.register_data['ghost'];sc=km.array([km.q(i) for i in range(ci.num_qubits) if km.q(i) not in set(x)]);sim.write_across_shots(x,bits(vals));sim.do(cf)
 assert np.array_equal(sim.read_across_shots(x),bits(want));assert not sim.read_across_shots(sc).any();mask=sim.read_across_shots(ghost);phase=(mask&bits(vals)).sum(0)%2!=0;assert np.array_equal(sim.read_phase_flipped_across_shots(),phase)
 sim.do(ci);assert np.array_equal(sim.read_across_shots(x),bits(vals));assert not sim.read_across_shots(sc).any();assert not sim.read_across_shots(ghost).any();assert not sim.read_phase_flipped_across_shots().any()
 for rev,c in [(False,cf),(True,ci)]:c.to_file(f'candidates/spooky_extension_inverse_compact24_{"reverse" if rev else "forward"}.kmb',format='kmb')
 r=dict(status='PASS',shots=shots,forward=dict(Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like),reverse=dict(Q=ci.num_qubits,T=ci.operation_statistics().num_toffoli_like),scope='in-place inverse pair; retained X-measurement phase checked before repair; zero, axes, random elements, clean quantum and classical scratch after pair')
 Path('results/spooky_extension_inverse_compact24_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()
