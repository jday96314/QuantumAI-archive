"""Repair one zero scalar pivot; a second zero in its pair stays unsupported."""
from functools import lru_cache
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder

@lru_cache(None)
def build(first=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());g=b.create_quantum_register(12,'g');cur=b.create_quantum_register(12,'cur');b.create_quantum_register(12,'next');prev=b.create_quantum_register(12,'prev_b');occupied=b.create_quantum_register(1,'occupied');out=b.create_quantum_register(1,'out');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 with m.opened(E,'current_pivot',12) as value:
  rt=m.prod_rt(value,[q.id for q in g],[q.id for q in cur])
  if not first:rt +=[('L',(q.id,),0,(r,)) for q,r in zip(prev,value)]
  E.mac(rt)
  with m.eq_test(E,value,0) as z:
   if first:E.cx(z,out[0].id)
   else:b.x(occupied);E.b.ccx(z,occupied[0].id,out[0].id);b.x(occupied)
  E.free(rt)
 ab.release_recycle('single_pivot_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
