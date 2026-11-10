"""Shared long-tail three-slot transform with first/second gap support.

Retain the high three raw moments, use the six-word cofactor frame, and
prepare P by one short middle product. The last quotient degree chooses
the existing shared shift-two/shift-three inverse. Every input restores;
decoding gaps crossing the last slot remain a separate scheduling task.
"""
from contextlib import ExitStack
import triple_gap_head as head
import triple_gap_frame as frame
from jacobi_decoder import cf
from variable_shift_inverse import apply as reciprocal

def block(m,E,math,W,inverse=False):
 H,R=W[:6],W[6:]
 if not R:head.apply(m,E,H,inverse);return
 with m.opened(E,'triple_raw_high',36) as rawbank:
  raw=[rawbank[12*i:12*(i+1)] for i in range(3)]
  def copy():
   for src,dst in zip(H[3:],raw):
    for a,b in zip(src,dst):E.cx(a,b)
  if inverse:head.apply(m,E,H,True)
  copy();head.apply(m,E,H)
  with m.opened(E,'triple_gap_charts',2) as flags:
   def census(back=False):
    if not back:
     with m.eq_test(E,H[0],0) as zero:E.cx(zero,flags[0])
    E.x(flags[0])
    with m.eq_test(E,H[2],0) as zero:E.b.ccx(flags[0],zero,flags[1])
    E.x(flags[0])
    if back:
     with m.eq_test(E,H[0],0) as zero:E.cx(zero,flags[0])
   census()
   with frame.parameters(m,E,H,flags) as (A,D,p):
    with ExitStack() as stack:
     P=[stack.enter_context(m.opened(E,f'triple_affine_{i}',12)) for i in range(3)];targets=[[] for _ in range(5)]
     for i,out in enumerate(P):targets[2+i]=[out]
     pr=m.polyprod_rt(A[1:],raw,targets);E.mac(pr)
     def shift():
      for src,dst in zip(P,R):
       for a,b in zip(src,dst):E.cx(a,b)
     if inverse:
      math.triangular(A,1,R,True);math.triangular(D,1,R);reciprocal(m,E,math,R,flags[1]);cf.scale(m,E,R,p)
      math.triangular(D,1,R,True);shift();math.triangular(A,1,R,True)
     else:
      math.triangular(A,1,R);shift();math.triangular(D,1,R);cf.scale(m,E,R,p,True);reciprocal(m,E,math,R,flags[1])
      math.triangular(D,1,R,True);math.triangular(A,1,R)
     E.free(pr)
   census(True)
  head.apply(m,E,H,True);copy()
  if not inverse:head.apply(m,E,H)
