"""Pad an even degree by x and normalize the explicit denominator lead."""
import kickmix as km
from extension_inverse import setup,Adapter
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create
from cubic_low_workspace_inverse import bind


def build(n,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());N=[b.create_quantum_register(12,f'N{i}') for i in range(n+2)];B=[b.create_quantum_register(12,f'B{i}') for i in range(n+1)];flag=b.create_quantum_register(1,'padding')[0];ab=Adapter(b);sc,E=lean_create(m,ab)
 def detect():b.flip_if_equal(B[-1],0,target=flag)
 def shift(back=False):
  for bank in [N,B]:
   for i in (range(len(bank)-1) if back else reversed(range(len(bank)-1))):b.cswap(flag,bank[i][0],bank[i+1][0]) if False else [b.cswap(flag,q,r) for q,r in zip(bank[i],bank[i+1])]
 if not inverse:detect();shift()
 lead=[q.id for q in B[-1]]
 with m.inv_word(E,lead) as li,m.opened(E,'raw_normalize_product',12) as tmp:
  for w in N[:n]+B[:n]:m.mon_slot(E,[q.id for q in w],lead if inverse else li,li if inverse else lead,tmp)
 if inverse:shift(True);detect()
 ab.release_recycle('raw_normalize_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
