"""Large prefix block using reversed-quotient cofactor coordinates.

Only upper raw moments are retained. The affine boundary reuses that bank
in place, arming a zero leading diagonal. Completion/singular decoding is
not asserted for this prefix policy; the map always restores under inversion.
"""
import lookahead_jacobi as look
import compact_lookahead_parameters as compact
import quotient_cofactor_frame as cofactor
from prefix_high_copy_jacobi import boundary
from variable_shift_inverse import apply as reciprocal
from jacobi_decoder import cf

def head(m,E,math,H,inverse=False):
 original=look.parameters;look.parameters=compact.parameters
 try:look.extract(m,E,math,H,inverse)
 finally:look.parameters=original

def block(m,E,math,W,pairs=6,inverse=False):
 assert pairs%2==0;H,R=W[:2*pairs],W[2*pairs:]
 if not R:head(m,E,math,H,inverse);return
 with m.opened(E,'prefix_quotient_raw_high',12*pairs) as rawbank:
  P=[rawbank[12*i:12*(i+1)] for i in range(pairs)]
  def copy():
   for src,dst in zip(H[pairs:],P):
    for a,b in zip(src,dst):E.cx(a,b)
  if inverse:head(m,E,math,H,True)
  copy();head(m,E,math,H)
  with m.opened(E,'prefix_quotient_last_gap',1) as last:
   def record():
    with m.eq_test(E,H[-4],0) as zero:E.cx(zero,last[0])
   record()
   with cofactor.parameters(m,E,H) as (A,D,p):
    boundary(m,E,A,P)
    def shift():
     for src,dst in zip(P,R):
      for a,b in zip(src,dst):E.cx(a,b)
    if inverse:
     math.triangular(A,1,R,True);math.triangular(D,1,R);reciprocal(m,E,math,R,last[0]);cf.scale(m,E,R,p)
     math.triangular(D,1,R,True);shift();math.triangular(A,1,R,True)
    else:
     math.triangular(A,1,R);shift();math.triangular(D,1,R);cf.scale(m,E,R,p,True);reciprocal(m,E,math,R,last[0])
     math.triangular(D,1,R,True);math.triangular(A,1,R)
    boundary(m,E,A,P,True)
   record()
  head(m,E,math,H,True);copy()
  if not inverse:head(m,E,math,H)
