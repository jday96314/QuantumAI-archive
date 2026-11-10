"""Blocked Jacobi transform with only the upper raw half-head retained.

Copy b raw moments, convert the 2b head in place to denominator continuants,
and form the boundary with one short polynomial product. The head is briefly
restored to erase the copy, then returned to Jacobi coordinates. This replaces
large symbolic affine formulas with extra short head transforms. Long-tail
arithmetic is shared. Interior-gap decoding is not claimed by this module.
"""
from contextlib import ExitStack
from jacobi_decoder import extract as head_extract,cf,pair
from continuant_head_frame import frame


def block(m,E,math,W,pairs=3,inverse=False):
 H,R=W[:2*pairs],W[2*pairs:]
 if not R:
  head_extract(m,E,math,H,inverse);return
 with m.opened(E,'raw_half_head',12*pairs) as bank:
  raw=[bank[12*i:12*(i+1)] for i in range(pairs)]
  def copy():
   for src,dst in zip(H[pairs:],raw):
    for a,b in zip(src,dst):E.cx(a,b)
  if inverse:head_extract(m,E,math,H,True)
  copy();head_extract(m,E,math,H)
  with m.opened(E,'high_copy_zero_flags',pairs-1) as flags:
   A,D,p=frame(m,E,H,flags)
   with ExitStack() as stack:
    P=[stack.enter_context(m.opened(E,f'high_copy_affine_{i}',12)) for i in range(pairs)]
    targets=[[] for _ in range(2*pairs-1)]
    for i,out in enumerate(P):targets[pairs-1+i]=[out]
    pr=m.polyprod_rt(A[1:],raw,targets);E.mac(pr)
    def shift():
     for src,dst in zip(P,R):
      for a,b in zip(src,dst):E.cx(a,b)
    def reciprocal():
     present=() in math.known_zero;math.known_zero.add(());math.inverse([(),(),*R])
     if not present:math.known_zero.remove(())
    if inverse:
     math.triangular(A,1,R,True);math.triangular(D,1,R);reciprocal();cf.scale(m,E,R,p)
     math.triangular(D,1,R,True);shift();math.triangular(A,1,R,True)
    else:
     math.triangular(A,1,R);shift();math.triangular(D,1,R);cf.scale(m,E,R,p,True);reciprocal()
     math.triangular(D,1,R,True);math.triangular(A,1,R)
    E.free(pr)
   frame(m,E,H,flags,True)
  head_extract(m,E,math,H,True);copy()
  if not inverse:head_extract(m,E,math,H)


def extract(m,E,math,W,inverse=False,pairs=3):
 for i in (reversed(range(0,len(W),2*pairs)) if inverse else range(0,len(W),2*pairs)):
  count=min(pairs,(len(W)-i)//2)
  if count==1:pair(m,E,math,W[i:],inverse)
  else:block(m,E,math,W[i:],count,inverse)
