"""In-place quadratic inversion with a held norm measurement.

The two coordinates scale in place. The old norm, rather than the old
24-bit source, carries the sole pending X-measurement phase.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup,extension_norm_rt,inverse_output_rt
from low_workspace_inverse import norm_rt as norm12_rt
from inverse_mask_kernel import phase_inverse
from tower_inverse import apply_map,map_ops


def init_delayed(m,E,src):
 names=[('ext_inv_N6',6),('ext_inv_M3',3),('ext_inv_S3',3),('ext_inv_Y6',6)];wins=[E.A.open(name,n) for name,n in names];n,mm,s,y=wins
 rN,rM,rS=m.norm_chain_terms(src,m.spec_rows(m.A0ROWS,src),m.spec_rows(m.A1ROWS,src),n,mm,s);rY=m.chain_rts(src,wins,src)[-2];E.macs([rN,rM,rS,rY]);E.frees([rN,rM,rS])
 for (name,_),word in reversed(list(zip(names[:3],wins[:3]))):E.A.close(name,word)
 out=E.A.open('ext_delayed_inverse',12);E.mac(inverse_output_rt(m,src,y,out));E.free(norm12_rt(m,out,y));E.A.close(names[3][0],y);return out


def inverse_mask_native(b,source,mask,fw,bk):
 word=km.array([km.q(q) for q in source]);bops=map_ops(bk);trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans)
 apply_map(b,word,bops);apply_map(b,mask,mops);phase_inverse(b,word,mask);apply_map(b,mask,mops,True);apply_map(b,word,bops,True)


def build(reverse=False):
 m,k,fw,bk,tr=setup();b=FastBuilder(km.CircuitBuilder());x=b.create_quantum_register(24,'x');ghost=b.create_classical_register(24,'ghost');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));src=[q.id for q in x]
 nn=E.A.open('extension_norm',12);E.mac(extension_norm_rt(m,src,nn,k));ni=init_delayed(m,E,nn)
 if reverse:
  for q,mask in zip(ni,ghost[:12]):ab.z(q,mask)
  b.bit_store(ghost,False)
 else:b.bit_store(ghost,False)
 b.cx(x[12:],x[:12])
 with m.opened(E,'extension_coordinate_scale',12) as temp:
  for word in [src[:12],src[12:]]:m.mon_slot(E,word,ni,nn,temp)
 if reverse:
  mask=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in nn]),mask);E.A.close('extension_norm',nn);inverse_mask_native(b,ni,mask,fw,bk);b.bit_store(mask,False);b.free(mask)
 else:
  b.hmr(km.array([km.q(q) for q in nn]),ghost[:12]);E.A.close('extension_norm',nn)
 E.free(extension_norm_rt(m,src,ni,k));E.A.close('ext_delayed_inverse',ni);ab.release_recycle('norm_inverse_emitter',sc);c=b.finish_circuit();held={v.id for v in ghost};return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits) if i not in held))


def test():
 cf=build();ci=build(True);m,k,fw,bk,tr=setup();rng=np.random.default_rng(142093);shots=16384;vals=rng.integers(0,1<<24,shots);vals[:4096]=np.arange(4096);vals[4096:8192]=np.arange(4096)<<12
 def norm(v):
  a,b=int(v)&4095,int(v)>>12;return m.gmul(a,a)^m.gmul(a,b)^m.gmul(k,m.gmul(b,b))
 def inv(v):
  a,b=int(v)&4095,int(v)>>12;ni=m.ginv(norm(v));return m.gmul(a^b,ni)|(m.gmul(b,ni)<<12)
 want=np.array([inv(v) for v in vals]);norms=np.array([norm(v) for v in vals]);bits=lambda v,w=24:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);rows=[]
 for seed in [142093,142094]:
  sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(ci);x=cf.register_data['x'];ghost=cf.register_data['ghost'];scratch=km.array([km.q(i) for i in range(max(cf.num_qubits,ci.num_qubits)) if km.q(i) not in set(x)]);sim.write_across_shots(x,bits(vals));sim.do(cf);mask=sim.read_across_shots(ghost);expected=((mask[:12]&bits(norms,12)).sum(0)%2).astype(bool)
  r=dict(seed=seed,value=bool((sim.read_across_shots(x)!=bits(want)).any()),scratch=bool(sim.read_across_shots(scratch).any()),phase=bool((sim.read_phase_flipped_across_shots()!=expected).any()));print(json.dumps(r),flush=True);rows.append(r);assert not any(r[k] for k in ['value','scratch','phase']);sim.do(ci);assert np.array_equal(sim.read_across_shots(x),bits(vals));assert not sim.read_across_shots(scratch).any();assert not sim.read_across_shots(ghost).any();assert not sim.read_phase_flipped_across_shots().any()
 sim=km.Simulator(batch_size=shots,seed=142095);sim.use_same_registers_as(ci);sim.write_across_shots(ci.register_data['x'],bits(vals));sim.do(ci);assert np.array_equal(sim.read_across_shots(ci.register_data['x']),bits(want));assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 for rev,c in [(False,cf),(True,ci)]:c.to_file(f'candidates/spooky_extension_norm_inverse_{"reverse" if rev else "forward"}.kmb',format='kmb')
 result=dict(status='PASS',shots=shots*3,forward=dict(Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like),reverse=dict(Q=ci.num_qubits,T=ci.operation_statistics().num_toffoli_like),rows=rows,scope='all two subfield axes and random extension inputs, zero, two forward/pair measurements plus fresh inverse, exact held norm phase, clean source/scratch/classical outcomes');Path('results/spooky_extension_norm_inverse_native.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':test()
