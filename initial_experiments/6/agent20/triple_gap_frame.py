"""Six-word cofactor frame for three linear slots or either degree-two gap.

Heads are normal(a,b,c,d,e,f), first gap(0,u,h,e,a,b), or second
gap(a,b,0,u,h,e). All return the same A/D/p wire layout. Second gap
has D's degree-two coefficient zero and uses a shift-three tail map.
This is only a parameter frame, not a complete three-slot extractor.
"""
from contextlib import contextmanager
import continuant_head_frame as regular
from quantum_control import selected
from jacobi_decoder import cf

def first(m,E,H,inverse=False):
 z,u,h,e,a,b=H
 def scale(back=False):
  with cf.armed(m,E,a),m.inv_word(E,a) as ai,m.opened(E,'triple_first_product',12) as t:
   m.mon_slot(E,u,ai if back else a,a if back else ai,t)
 def route(back=False):
  for x,y in ([(u,a),(z,u)] if back else [(z,u),(u,a)]):E.swaps(None,list(zip(x,y)))
 def change_b():
  for x,y in zip(e,b):E.cx(x,y)
 def change_z():E.xor(m.prod_rt(z,b,e)+[('L',(x,),0,(y,)) for x,y in zip(h,z)])
 def change_a():E.xor(m.prod_rt(a,b,h))
 if inverse:route(True);change_b();change_z();change_a();scale(True)
 else:scale();change_a();change_z();change_b();route()

def second(m,E,H,inverse=False):
 a,b,z,u,h,e=H
 def scale(back=False):
  with cf.armed(m,E,u),m.inv_word(E,u) as ui,m.opened(E,'triple_second_product',12) as t:
   m.mon_slot(E,a,ui if back else u,u if back else ui,t)
 def change_u():E.xor(m.prod_rt(u,b,h))
 def change_h():E.xor(m.prod_rt(h,b,e))
 def change_e():
  for x,y in zip(b,e):E.cx(x,y)
 if inverse:
  E.swaps(None,list(zip(b,u)));change_e();change_h();change_u();scale(True)
 else:scale();change_u();change_h();change_e();E.swaps(None,list(zip(b,u)))

@contextmanager
def parameters(m,E,H,flags):
 first_gap,second_gap=flags
 with m.opened(E,'triple_normal_completion',2) as completion:
  def normal(back=False):
   with m.eq_test(E,flags,0) as active,selected(m,E,active) as (cm,ce):regular.frame(cm,ce,H,completion,back)
  def branch(fn,flag,back=False):
   with selected(m,E,flag) as (cm,ce):fn(cm,ce,H,back)
  normal();branch(first,first_gap);branch(second,second_gap)
  _,A,D,p=regular.layout(H);yield A,D,p
  branch(second,second_gap,True);branch(first,first_gap,True);normal(True)
