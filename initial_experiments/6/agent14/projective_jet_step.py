"""Reversible Mobius action on a two-term projective infinity jet.

The infinity chart stores the reciprocal jet with zero constant coefficient.
This handles an isolated scalar pole without a degree counter or jet gauge.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create
from schur_control import caxpy
from spooky_extension_norm_inverse import init_delayed
from cubic_mask_phase_joint import native_phase_inverse
from tower_inverse import apply_map,map_ops


def build(inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());a=b.create_quantum_register(12,'a');bb=b.create_quantum_register(12,'b');h=b.create_quantum_register(12,'h');k=b.create_quantum_register(12,'k');tag=b.create_quantum_register(1,'infinity')[0];ab=Adapter(b);sc,E=lean_create(m,ab);H=[q.id for q in h];K=[q.id for q in k];z=b.alloc_qubits(1)[0];sq=map_ops([m.gmul(1<<i,1<<i) for i in range(12)])
 def translate():
  b.x(tag);caxpy(m,E,tag.id,[q.id for q in bb],H,1);caxpy(m,E,tag.id,[q.id for q in a],K,1);b.x(tag)
 def detect():
  with m.eq_test(E,H,0) as zz:ab.cx(zz,z.id)
 def reciprocal():
  out=init_delayed(m,E,H);reg=km.array([km.q(q) for q in out]);apply_map(b,h,sq);apply_map(b,reg,sq)
  with m.opened(E,'projective_jet_scale',12) as tmp:m.mon_slot(E,K,out,H,tmp)
  apply_map(b,reg,sq,True);apply_map(b,h,sq,True);b.cx(reg,k);mask=b.alloc_bits(12);b.hmr(h,mask)
  for q,r in zip(H,out):ab.swap(q,r)
  E.A.close('ext_delayed_inverse',out);native_phase_inverse(b,H,mask);b.bit_store(mask,False);b.free(mask)
 if inverse:
  detect();b.cx(z,tag);b.cx(z,h[0]);b.cx(z,k[0]);reciprocal();b.cx(z,h[0]);detect();translate()
 else:
  translate();detect();b.cx(z,h[0]);reciprocal();b.cx(z,k[0]);b.cx(z,h[0]);b.cx(z,tag);detect()
 b.free(km.array([z]));ab.release_recycle('projective_jet_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test():
 m=setup()[0];rng=np.random.default_rng(142118);shots=8192;aa,bb,hh,kk=rng.integers(0,4096,(4,shots));tag=np.repeat([0,1],4096);hh[:4096]=np.arange(4096);hh[4096:]=0;bb[::4]=hh[::4];wantH=[];wantK=[];wantT=[]
 for a,b,h,k,t in zip(aa,bb,hh,kk,tag):
  if t:outH=0;outK=int(k);outT=0
  else:
   d=int(h)^int(b);e=int(k)^int(a);outT=int(not d);outH=m.ginv(d);outK=e if not d else m.gmul(e,m.gmul(outH,outH))^outH
  wantH.append(outH);wantK.append(outK);wantT.append(outT)
 bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool);cf=build();ci=build(True);src={'a':aa,'b':bb,'h':hh,'k':kk,'infinity':tag};dst={'a':aa,'b':bb,'h':wantH,'k':wantK,'infinity':wantT};rows=[]
 for inv,c,s,d in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142118+inv);sim.use_same_registers_as(c)
  for key,val in s.items():sim.write_across_shots(c.register_data[key],bits(val,len(c.register_data[key])))
  sim.do(c);errors=[key for key,val in d.items() if not np.array_equal(sim.read_across_shots(c.register_data[key]),bits(val,len(c.register_data[key])))];pub={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in pub]);r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
  sim.do(cf if inv else ci)
  for key,val in s.items():assert np.array_equal(sim.read_across_shots(c.register_data[key]),bits(val,len(c.register_data[key])))
  assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any()
 Path('results/projective_jet_step_native.json').write_text(json.dumps(dict(status='PASS',shots=shots*2,rows=rows),indent=2)+'\n')
if __name__=='__main__':test()
