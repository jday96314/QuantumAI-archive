"""Narrow early blocks using an invertible upper-triangular boundary map.

The retained high half-head is changed in place into its affine boundary.
The diagonal and scalar heads are armed at zero. Thus the circuit always
restores under inversion; decoding is asserted only for nonzero early
Jacobi scalar heads and a nonzero block-denominator leading coefficient.
A fixed early-block schedule confines completion exceptions to low weights.
Later blocks use the full completion-aware high-copy implementation.
"""
from jacobi_decoder import extract as head_extract,cf,pair
from continuant_head_frame import layout
import high_copy_jacobi as late


def simple_frame(m,E,H,inverse=False):
 frames,A,D,p=layout(H)
 for P,B,beta,alpha in (reversed(frames) if inverse else frames):
  def high(back=False):
   if back:
    for a,b in zip(B[-1],alpha):E.cx(a,b)
    if len(B)>1:
     for a,b in zip(B[-2],beta):E.cx(a,b)
    E.xor(m.prod_rt(beta,alpha,B[-1]))
   else:
    E.xor(m.prod_rt(beta,alpha,B[-1]))
    if len(B)>1:
     for a,b in zip(B[-2],beta):E.cx(a,b)
    for a,b in zip(B[-1],alpha):E.cx(a,b)
  if inverse:high(True)
  with cf.armed(m,E,beta),m.inv_word(E,beta) as bi,m.opened(E,'prefix_frame_product',12) as t:
   if not inverse:m.mon_slot(E,p,beta,bi,t)
   for i,w in enumerate(P):
    if not inverse:m.mon_slot(E,w,beta,bi,t)
    E.xor(m.prod_rt(w,alpha,B[i]))
    if i:
     for a,b in zip(B[i-1],w):E.cx(a,b)
    if inverse:m.mon_slot(E,w,bi,beta,t)
   if inverse:m.mon_slot(E,p,bi,beta,t)
  if not inverse:high()
 return A,D,p


def boundary(m,E,A,P,inverse=False):
 b=len(P)
 with cf.armed(m,E,A[-1]),m.inv_word(E,A[-1]) as ai,m.opened(E,'prefix_boundary_product',12) as tmp:
  for i in (reversed(range(b)) if inverse else range(b)):
   if not inverse:m.mon_slot(E,P[i],A[-1],ai,tmp)
   for k in range(i+1,b):E.xor(m.prod_rt(P[i],A[b+i-k],P[k]))
   if inverse:m.mon_slot(E,P[i],ai,A[-1],tmp)


def block(m,E,math,W,pairs=6,inverse=False):
 H,R=W[:2*pairs],W[2*pairs:]
 if not R:head_extract(m,E,math,H,inverse);return
 with m.opened(E,'prefix_high_head',12*pairs) as bank:
  P=[bank[12*i:12*(i+1)] for i in range(pairs)]
  def copy():
   for src,dst in zip(H[pairs:],P):
    for a,b in zip(src,dst):E.cx(a,b)
  if inverse:head_extract(m,E,math,H,True)
  copy();head_extract(m,E,math,H);A,D,p=simple_frame(m,E,H);boundary(m,E,A,P)
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
  boundary(m,E,A,P,True);simple_frame(m,E,H,True);head_extract(m,E,math,H,True);copy()
  if not inverse:head_extract(m,E,math,H)


def schedule(n,prefix=6,through=18,suffix=2):
 counts=[];k=0
 while k<n:
  b=min(prefix,through-k,n-k) if k<through else min(suffix,n-k)
  counts.append(b);k+=b
 return counts


def extract(m,E,math,W,inverse=False,counts=None,through=18):
 counts=counts or schedule(len(W)//2,through=through);stages=[];k=0
 for b in counts:stages.append((k,b));k+=b
 assert k==len(W)//2
 for k,b in (reversed(stages) if inverse else stages):
  if b==1:pair(m,E,math,W[2*k:],inverse)
  elif k<through:block(m,E,math,W[2*k:],b,inverse)
  else:late.block(m,E,math,W[2*k:],b,inverse)
