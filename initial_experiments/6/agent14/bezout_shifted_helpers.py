"""Monomial-gcd completion and x+1 companion repair helpers."""
from contextlib import contextmanager
import kickmix as km


def translate(m,E,words,monic=False):
 degree=len(words) if monic else len(words)-1
 for i in range(len(words)):
  for j in range(i+1,len(words)):
   if i&j==i:
    for q,r in zip(words[j],words[i]):E.cx(q,r)
  if monic and i&degree==i:E.x(words[i][0])


@contextmanager
def arm(m,E,g,degree,history,padding):
 index,flag=history;b=E.b.raw
 def toggle():
  b.flip_if_less_than(km.array([km.q(q) for q in padding]),2*degree+1,target=km.q(flag));E.x(flag);
  for j,word in enumerate(g):
   if j&(len(g)-1)==j:E.cx(flag,word[0])
  E.x(flag);b.flip_if_less_than(km.array([km.q(q) for q in padding]),2*degree+1,target=km.q(flag))
 toggle()
 try:yield
 finally:toggle()
