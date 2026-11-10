"""One-product lookahead, valid for a nonzero regular or reciprocal frame."""
from functools import lru_cache
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from schur_control_templates import finish


@lru_cache(None)
def flagger(first=False,reciprocal=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());g=b.create_quantum_register(12,'g');cur=b.create_quantum_register(12,'cur');nxt=b.create_quantum_register(12,'next');pb=b.create_quantum_register(12,'prev_b');prev=b.create_quantum_register(1,'prev_flag');out=b.create_quantum_register(1,'out_flag');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));G=[q.id for q in g];C=[q.id for q in cur];P=[q.id for q in pb]
 with m.opened(E,'current_zero_test',12) as bc:
  rt=[('L',(q,),0,(r,)) for q,r in zip(C,bc)] if first or reciprocal else m.prod_rt(bc,G,C)
  if not first:
   rt+=m.prod_rt(bc,G,P) if reciprocal else [('L',(q,),0,(r,)) for q,r in zip(P,bc)]
  E.mac(rt);b.cx(cur,nxt)
  with m.eq_test(E,bc,0) as zc,m.eq_test(E,[q.id for q in nxt],0) as zn:
   E.x(zn)
   if first:E.b.ccx(zc,zn,out[0].id)
   else:
    with m.opened(E,'lookahead_and',1) as temp:
     rr=[('P',(zc,),0,(zn,),0,(temp[0],))];E.mac(rr);b.x(prev);E.b.ccx(temp[0],prev[0].id,out[0].id);b.x(prev);E.free(rr)
   E.x(zn)
  b.cx(cur,nxt);E.free(rt)
 return finish(b,ab,sc)
