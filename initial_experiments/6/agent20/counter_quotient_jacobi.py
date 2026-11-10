"""Completion-aware quotient synthesis using log(n)+1 metadata bits.

Count zero Jacobi scalar heads without retaining n predicates. Arm only
index d=n-count, and retain one bit saying whether that scalar was zero.
The count and armed bit suffice for literal inversion on arbitrary heads;
on regular/complete inputs they select exactly the true prefix boundary.
The controlled increment uses k-1 ANDs and phase-correct measurement erasure
before the corresponding low source bits are changed.
"""
from contextlib import contextmanager
import quotient_jacobi_synthesis as qj
import completed_quotient_jacobi as prior


@contextmanager
def equal(m,E,bits,value):
 if len(bits)==1:
  if not value:E.x(bits[0])
  yield bits[0]
  if not value:E.x(bits[0])
 else:
  with m.eq_test(E,bits,value) as match:yield match


def increment(m,E,bits,control,inverse=False):
 if inverse:
  for q in bits:E.x(q)
 with m.opened(E,'quotient_degree_increment',max(0,len(bits)-1)) as carry:
  terms=[]
  for j,target in enumerate(carry):
   rt=m.one_rt((control if j==0 else carry[j-1],),0,(bits[j],),0,target);terms.append(rt);E.mac(rt)
  for j in reversed(range(1,len(bits))):E.cx(carry[j-1],bits[j]);E.free(terms[j-1])
  E.cx(control,bits[0])
 if inverse:
  for q in bits:E.x(q)


def synthesize(m,E,W,c,metadata,inverse=False):
 n=len(W)//2;count,armed=metadata[:-1],metadata[-1];F,G=qj.layout(W,c)
 def census(back=False):
  for a in (reversed(W[::2]) if back else W[::2]):
   with m.eq_test(E,a,0) as zero:increment(m,E,count,zero,back)
 def arm_record():
  for i,a in enumerate(W[::2]):
   with equal(m,E,count,n-i) as match,m.eq_test(E,a,0) as zero:E.b.ccx(match,zero,armed)
 def flip():
  for i,a in enumerate(W[::2]):
   with equal(m,E,count,n-i) as match:E.b.ccx(match,armed,a[0])
 def terminal():
  with equal(m,E,count,0) as full:E.cx(full,c[0])
 if inverse:
  qj.poly.extract(m,E,F,G);flip();arm_record();terminal();qj.convert(m,E,W,True);census(True)
 else:
  census();qj.convert(m,E,W);terminal();arm_record();flip();qj.poly.extract(m,E,F,G,True)
 return F,G


def model(m,words):
 F,G,flags=prior.model(m,words);n=len(words)//2;zeros=sum(flags)
 return F,G,[(zeros>>k)&1 for k in range(n.bit_length())]+[int(zeros>0)]
