"""Consume the final remainder directly into the last Euclid cofactors.

When the remainder has exactly2b words, there is no tail to update. Its words
can hold the quotient transform directly; the duplicate head and the two
speculative cofactor passes of the general block are unnecessary.
"""
from blocked_monic_euclid import schedule,cofactor
from mobius_series import apply as mobius
from continued_fraction_series import apply as cf
from deferred_euclid_padding import fake_undo_ops,common_padding_ops
from controlled_linear import apply as controlled_linear


def apply(m,E,math,Bp,B,C,pinv,width,inverse=False,flags=None):
 assert len(C)==2*width and len(flags)==width
 Q=list(reversed(C));ff,(fp,fb)=schedule(Bp,B,Q);rev=[B[0],*reversed(B)];drev=[B[0],*reversed(Bp)]
 def head(inv=False):
  if inv:cf(m,E,math,Q,True,flags,True)
  mobius(m,E,math,rev,drev,pinv,Q,inv,p_is_inverse=True)
  if not inv:cf(m,E,math,Q,False,flags,True)
 def cof(inv=False):
  for j in (reversed(range(width)) if inv else range(width)):
   bp,bb,beta,alpha=ff[j];ops=fake_undo_ops(bp,bb,beta,alpha)+common_padding_ops(bp,bb,[beta,alpha],False)
   if inv:controlled_linear(E.b,ops,flags[j],True)
   cofactor(m,E,bp,bb,beta,alpha,pinv,inv,True)
   if not inv:controlled_linear(E.b,ops,flags[j])
 if inverse:cof(True);head(True)
 else:head();cof()
 return fp,fb,[],pinv
