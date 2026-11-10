"""Reusable native templates for controlled node operations and lookahead."""
from functools import lru_cache
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup
from schur_control import caxpy,cimul


def finish(b,ab,sc):
 ab.release_recycle('control_template_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))

@lru_cache(None)
def axpy(c):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());flag=b.create_quantum_register(1,'flag');src=b.create_quantum_register(12,'src');dst=b.create_quantum_register(12,'dst');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));caxpy(m,E,flag[0].id,[q.id for q in src],[q.id for q in dst],c);return finish(b,ab,sc)

@lru_cache(None)
def scale(c,inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());flag=b.create_quantum_register(1,'flag');src=b.create_quantum_register(12,'src');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));cimul(m,E,flag[0].id,[q.id for q in src],c,inverse);return finish(b,ab,sc)

@lru_cache(None)
def flagger(first=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());g=b.create_quantum_register(12,'g');cur=b.create_quantum_register(12,'cur');nxt=b.create_quantum_register(12,'next');pb=b.create_quantum_register(12,'prev_b');prev=b.create_quantum_register(1,'prev_flag');out=b.create_quantum_register(1,'out_flag');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));G=[q.id for q in g];C=[q.id for q in cur];N=[q.id for q in nxt];P=[q.id for q in pb]
 with m.opened(E,'current_actual_b',12) as bc,m.opened(E,'next_actual_b',12) as bn:
  rc=m.prod_rt(bc,G,C);rn=m.prod_rt(bn,G,N)
  if not first:
   rc+=[('L',(q,),0,(v,)) for q,v in zip(P,bc)];rn+=[('L',(q,),0,(v,)) for q,v in zip(P,bn)]
  E.mac(rc);E.mac(rn)
  with m.eq_test(E,bc,0) as zc,m.eq_test(E,bn,0) as zn:
   E.x(zn)
   if first:E.b.ccx(zc,zn,out[0].id)
   else:
    with m.opened(E,'lookahead_and',1) as temp:
     rr=[('P',(zc,),0,(zn,),0,(temp[0],))];E.mac(rr);b.x(prev);E.b.ccx(temp[0],prev[0].id,out[0].id);b.x(prev);E.free(rr)
   E.x(zn)
  E.free(rn);E.free(rc)
 return finish(b,ab,sc)
