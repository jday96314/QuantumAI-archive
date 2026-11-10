"""Reversible completion-aware conversion from Schur heads to affine factors.

Retains the completed-prefix length, a boundary reciprocal, and an endpoint
ratio. These outputs are explicit; their width and cleanup are charged.
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
 m=setup()[0];m.inv_word=bind(m);n=N//2;b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(N)];degree=b.create_quantum_register(n.bit_length(),'degree');seed=b.create_quantum_register(12,'seed');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 def count(back=False):
  for i in (reversed(range(n)) if back else range(n)):
   with m.eq_test(E,[q.id for q in S[2*i+1]],0) as z:
    E.x(z);add(b,km.array([km.q(z)]),degree,subtract=back);E.x(z)
 def armor(back=False):
  flag=b.alloc_qubits(1)[0]
  if not back:b.flip_if_equal(degree,0,target=flag);b.cx(flag,S[1][0]);b.flip_if_equal(degree,0,target=flag)
  for i in (reversed(range(1,n)) if back else range(1,n)):
   b.flip_if_less_than(degree,i+1,target=flag)
   for j in range(24):b.ccx(flag,S[2*i-2+j//12][j%12],S[2*i+j//12][j%12])
   b.flip_if_less_than(degree,i+1,target=flag)
  if back:b.flip_if_equal(degree,0,target=flag);b.cx(flag,S[1][0]);b.flip_if_equal(degree,0,target=flag)
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
 if inverse:factors(True);endpoint(True);armor(True);count(True)
 else:count();armor();endpoint();factors()
 ab.release_recycle('affine_parameter_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


