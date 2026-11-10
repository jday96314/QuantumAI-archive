"""Approximate completion chart with explicit X-outcome bookkeeping.

Coincident nonzero active heads and accidental identity affine factors are
unsupported rare charts. Complete native verification measures their cost.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from base_inverse import inplace as invert12
from short_adder import add
from schur_control import caxpy


def build(N,inverse=False):
 m=setup()[0];m.inv_word=bind(m);n=N//2;b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(N)];seed=b.create_quantum_register(12,'seed');ghost=b.create_classical_register(n,'completion_ghost');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 def armor(back=False):
  flag=b.alloc_qubits(1)[0]
  for i in (reversed(range(n)) if back else range(n)):
   current=km.array(list(S[2*i])+list(S[2*i+1]));prev=km.array(list(S[2*i-2])+list(S[2*i-1])) if i else None
   if back:
    if i:b.cx(prev,current)
    b.flip_if_equal(current,0 if i else 1<<12,target=flag)
    if i:b.cx(prev,current)
    b.push_condition(ghost[i]);b.z(flag);b.pop_condition();b.bit_store(ghost[i],False)
   else:b.flip_if_equal(current,0,target=flag)
   if i:
    for a,c in zip(prev,current):b.ccx(flag,a,c)
   else:b.cx(flag,S[1][0])
   if back:b.flip_if_equal(current,0,target=flag)
   else:b.bit_store(ghost[i],False);b.hmr(flag,ghost[i])
  b.free(km.array([flag]))
 def endpoint(back=False):
  with m.inv_word(E,[q.id for q in S[-1]]) as bi:
   rt=m.prod_rt([q.id for q in seed],[q.id for q in S[-2]],bi)
   if back:E.free(rt)
   else:E.mac(rt)
 def factors(back=False):
  for i in (range(n-1) if back else reversed(range(n-1))):
   olda,oldb=S[2*i:2*i+2];na,nb=S[2*i+2:2*i+4];B=[q.id for q in oldb]
   if not back:b.cx(olda,na)
   with m.inv_word(E,B) as bi,m.opened(E,'affine_parameter_product',12) as temp:
    for w in [na,nb]:m.mon_slot(E,[q.id for q in w],B if back else bi,bi if back else B,temp)
   if back:b.cx(olda,na)
 if inverse:factors(True);endpoint(True);armor(True)
 else:armor();endpoint();factors()
 ab.release_recycle('affine_parameter_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(n,c.num_bits)))


