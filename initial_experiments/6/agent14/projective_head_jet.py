"""Consume (a,b) into the infinity jet of q/(a+b*x), in place."""
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from lean_emitter import create
from spooky_extension_norm_inverse import init_delayed
from cubic_mask_phase_joint import native_phase_inverse
from tower_inverse import apply_map,map_ops


def build(inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());k=b.create_quantum_register(12,'k');h=b.create_quantum_register(12,'h');tag=b.create_quantum_register(1,'infinity')[0];ab=Adapter(b);sc,E=create(m,ab);H=[q.id for q in h];K=[q.id for q in k];sq=map_ops([m.gmul(1<<i,1<<i) for i in range(12)])
 def detect():
  with m.eq_test(E,H,0) as z:ab.cx(z,tag.id)
 def recip():
  out=init_delayed(m,E,H);reg=km.array([km.q(q) for q in out]);apply_map(b,h,sq);apply_map(b,reg,sq)
  with m.opened(E,'head_jet_scale',12) as tmp:m.mon_slot(E,K,out,H,tmp)
  apply_map(b,reg,sq,True);apply_map(b,h,sq,True);b.cx(reg,k);mask=b.alloc_bits(12);b.hmr(h,mask)
  for q,r in zip(H,out):ab.swap(q,r)
  E.A.close('ext_delayed_inverse',out);native_phase_inverse(b,H,mask);b.bit_store(mask,False);b.free(mask)
 if inverse:b.cx(tag,h[0]);b.cx(tag,k[0]);recip();b.cx(tag,h[0]);detect()
 else:detect();b.cx(tag,h[0]);recip();b.cx(tag,k[0]);b.cx(tag,h[0])
 ab.release_recycle('head_jet_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
