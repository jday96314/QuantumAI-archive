"""Completion chart inferred from two leading coefficients.

An actual zero Euclidean polynomial is armed by (x+1)^degree. The marker
can be erased against the two leading coefficients of the armed polynomial.
Accidental matches and gaps exceeding one are unsupported approximate charts.
The entrance/exit are literal inverse sequences, without a padding counter.
"""
from contextlib import contextmanager
from bezout_shifted_helpers import translate


@contextmanager
def arm(m,E,g,degree,history,padding=None,f=None):
 if len(g)<2:
  yield
  return
 index,flag=history;n=len(g)-1
 def before():
  with m.eq_test(E,g[-1]+g[-2],0) as z:E.cx(z,flag)
 def after():
  code=1 | ((n&1)<<12)
  assert f is not None
  if len(f)==degree+1:
   if degree&1:
    for q,r in zip(f[-1],f[-2]):E.cx(q,r)
   with m.eq_test(E,g[-1]+g[-2]+f[-2],code) as z:E.cx(z,flag)
   if degree&1:
    for q,r in zip(f[-1],f[-2]):E.cx(q,r)
  else:
   assert len(f)==degree
   with m.eq_test(E,g[-1]+g[-2]+f[-1],code|((degree&1)<<24)) as z:E.cx(z,flag)
 def toggle():
  for j,word in enumerate(g):
   if j&n==j:E.cx(flag,word[0])
 before();toggle();after()
 try:yield
 finally:after();toggle();before()
