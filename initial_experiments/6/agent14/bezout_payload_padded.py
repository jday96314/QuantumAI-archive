"""Regularized Bezout action on a consistent monomial-gcd chart."""
from contextlib import nullcontext
import kickmix as km
from bezout_payload_batched import divide
from bezout_payload_shared import companion as raw_companion
from bezout_padding_helpers import translate,arm


def companion(m,E,f,payloads,flag,inverse=False):
 for words in [f]+payloads:translate(m,E,words)
 raw_companion(m,E,f,payloads,flag,inverse)
 for words in [f]+payloads:translate(m,E,words)


def apply(m,E,f,g,payloads,inverse=False,g_one=False,history=None,padding=None):
 mm,n=len(f)-1,len(g)-1;assert all(len(v)==mm for v in payloads) and mm>=n>=0
 with arm(m,E,g,mm,history,padding) if mm>n else nullcontext():
  if mm==n:
   divide(m,E,f,g,one=g_one);apply(m,E,g,f[:n],payloads,inverse,history=history,padding=padding);divide(m,E,f,g,True,g_one);return
  assert mm==n+1
  if not n:
   with m.opened(E,'bez_terminal_pivot',12) as d:
    rt=m.prod_rt(d,f[-1],g[-1]);E.mac(rt)
    with m.inv_word(E,d) as di,m.opened(E,'bez_terminal_scale',12) as tmp:
     for y in payloads:m.mon_slot(E,y[-1],d if inverse else di,di if inverse else d,tmp)
    E.free(rt)
   return
  b=E.b.raw;index,flag=history;code=mm
  def detect():
   with m.eq_test(E,index,0) as z0,m.eq_test(E,g[-1],0) as zg:E.b.ccx(z0,zg,flag)
  def store():
   for i,q in enumerate(index):
    if code>>i&1:E.cx(flag,q)
  def query():b.flip_if_equal(km.array([km.q(q) for q in index]),code,target=km.q(flag))
  def rotate(back=False):
   if back:
    for i in reversed(range(n)):
     for q,r in zip(g[i+1],g[i]):b.ccx(km.q(flag),km.q(q),km.q(r))
   for i in (range(n) if back else reversed(range(n))):
    for q,r in zip(g[i],g[i+1]):b.cswap(km.q(flag),km.q(q),km.q(r))
   if not back:
    for i in range(n):
     for q,r in zip(g[i+1],g[i]):b.ccx(km.q(flag),km.q(q),km.q(r))
  detect();store()
  if not inverse:companion(m,E,f,payloads,flag)
  rotate()
  query()
  if True:
   def scale_last(scale,eraser):
    with m.opened(E,'bez_fast_scale',12) as tmp:
     for y in payloads:m.mon_slot(E,y[-1],scale,eraser,tmp)
   def cross(post=False):
    for i in range(n):
     for y in payloads:E.mac(m.prod_rt(y[-1] if post else y[i],g[i],y[i] if post else y[-1]))
   if not n:
    with m.opened(E,'bez_fast_pivot',12) as d:
     rt=m.prod_rt(d,f[-1],g[-1]);E.mac(rt)
     with m.inv_word(E,d) as di:scale_last(d if inverse else di,di if inverse else d)
     E.free(rt)
   else:
    with m.inv_word(E,g[-1]) as gi:
     scale_last(g[-1] if inverse else gi,gi if inverse else g[-1]);cross(inverse);divide(m,E,f,g,provided_inverse=gi)
    if inverse:apply(m,E,g,f[:n],[y[:n] for y in payloads],True,history=history,padding=padding)
    with m.inv_word(E,f[-1]) as qi:scale_last(f[-1] if inverse else qi,qi if inverse else f[-1])
    if not inverse:apply(m,E,g,f[:n],[y[:n] for y in payloads],history=history,padding=padding)
    with m.inv_word(E,g[-1]) as gi:
     divide(m,E,f,g,True,provided_inverse=gi);cross(not inverse);scale_last(g[-1] if inverse else gi,gi if inverse else g[-1])
  query();rotate(True)
  if inverse:companion(m,E,f,payloads,flag,True)
  store();detect()
