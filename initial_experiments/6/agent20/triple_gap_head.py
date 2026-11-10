"""Shared six-moment head covering either of the first two gap slots.

Route one common four-moment lookahead head. For the leading-gap chart,
only two residual moments remain; their Mobius update collapses to two
scalings and short affine corrections. This avoids nested controlled blocks.
Every zero pivot is armed and the inverse is literal on all inputs.
"""
import lookahead_jacobi as look
from routed_quantum_control import selected
from compact_blocked_cfraction import Arithmetic
from jacobi_decoder import cf

def scale(m,E,words,pivot,flag,inverse=False):
 with cf.armed(m,E,pivot),m.inv_word(E,pivot) as pi,m.opened(E,'triple_head_scale_product',12) as tmp,m.opened(E,'triple_head_scale_routed',12) as routed:
  for word in words:
   E.swaps(flag,list(zip(word,routed)))
   m.mon_slot(E,routed,pi if inverse else pivot,pivot if inverse else pi,tmp)
   E.swaps(flag,list(zip(word,routed)))

def pair(m,E,H,flag,inverse=False):
 def reciprocal():
  with selected(m,E,flag) as (cm,ce):Arithmetic(cm,ce).inverse(H)
 if inverse:reciprocal();scale(m,E,H[1:],H[0],flag)
 else:scale(m,E,H[1:],H[0],flag,True);reciprocal()

def gap_tail(m,E,R,H,flag,inverse=False):
 _,u,h,e=H;math=Arithmetic(m,E)
 with m.opened(E,'triple_head_gap_active',1) as active:
  def toggle():
   with m.eq_test(E,u,0) as zero:
    E.x(zero);E.b.ccx(flag,zero,active[0]);E.x(zero)
  def affine():
   with m.opened(E,'triple_head_affine',12) as out:
    cr=math.cube_rts(m.spec_plain(e),lambda bits:[out[j] for j in bits]);E.mac(cr)
    for a,b in zip(out,R[0]):E.b.ccx(active[0],a,b)
    E.free(cr)
    rt=[];m.bil_terms(rt,m.KT,m.spec_plain(h),math.spec([(1,e,2)]),tgt=out)
    m.lin_rows_terms(rt,math.rows(1,2),h,out);m.lin_rows_terms(rt,math.rows(1,4),e,out);E.mac(rt)
    for a,b in zip(out,R[1]):E.b.ccx(active[0],a,b)
    E.free(rt)
  toggle()
  if inverse:affine();scale(m,E,R,u,flag)
  else:scale(m,E,R,u,flag,True);affine()
  toggle()

def apply(m,E,H,inverse=False):
 assert len(H)==6
 with look.zero_head(m,E,H[0]) as first_gap:
  def normal_pair(back=False):
   E.x(first_gap);pair(m,E,H,first_gap,back);E.x(first_gap)
  def route(back=False):
   swaps=[(H[:2],H[4:]),(H[2:4],H[4:])]
   for A,B in (reversed(swaps) if back else swaps):E.swaps(first_gap,list(zip(sum(A,[]),sum(B,[]))))
  def common_head(back=False):
   with look.zero_head(m,E,H[2]) as flag:
    def normal_coordinates(undo=False):
     E.x(flag)
     with selected(m,E,flag) as (cm,ce):
      if not undo:
       for a,b in zip(H[3],H[5]):ce.cx(a,b)
      ce.xor(cm.prod_rt(H[4],H[3],H[5]))
      if undo:
       for a,b in zip(H[3],H[5]):ce.cx(a,b)
     E.x(flag)
    if back:normal_coordinates(True);look.head(m,E,H[2:],flag,True)
    else:look.head(m,E,H[2:],flag);normal_coordinates()
  if inverse:
   route();pair(m,E,H[:2],first_gap,True);gap_tail(m,E,H[:2],H[2:],first_gap,True);common_head(True);route(True);normal_pair(True)
  else:
   normal_pair();route();common_head();gap_tail(m,E,H[:2],H[2:],first_gap);pair(m,E,H[:2],first_gap);route(True)
