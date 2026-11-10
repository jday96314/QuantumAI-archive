"""Compact Euclid through a full local-ratio series instead of a copied head.

Convert the retained remainder into its ratio, take b Schur quotients, consume
those same words into the cofactors once, and convert the surviving ratio back
to a remainder. All field data stays in the original N-word state.
"""
from blocked_monic_euclid import schedule,cofactor
from mobius_series import apply as mobius
from deferred_euclid_padding import fake_undo_ops,common_padding_ops
from controlled_linear import apply as controlled_linear


def schur(m,E,math,F,count,flags,inverse=False):
 if not count:return
 if inverse:schur(m,E,math,F[2:],count-1,flags[1:],True)
 if not inverse:
  with m.eq_test(E,F[0],0) as zero:E.cx(zero,flags[0])
  E.cx(flags[0],F[0][0])
 with m.inv_word(E,F[0]) as iv,m.opened(E,'series_schur_product',12) as prod:
  if not inverse:
   for word in F[1:]:m.mon_slot(E,word,iv,F[0],prod)
  math.inverse(F)
  if inverse:
   for word in F[1:]:m.mon_slot(E,word,F[0],iv,prod)
 if inverse:
  E.cx(flags[0],F[0][0])
  with m.eq_test(E,F[0],0) as zero:E.cx(zero,flags[0])
 if not inverse:schur(m,E,math,F[2:],count-1,flags[1:])


def apply(m,E,math,Bp,B,C,pinv,width,inverse=False,flags=None):
 assert len(C)>=2*width and len(flags)==width
 F=list(reversed(C));Q=F[:2*width];tail=F[2*width:];ff,(fp,fb)=schedule(Bp,B,Q)
 oldrev=[B[0],*reversed(B)];olddrev=[B[0],*reversed(Bp)];newrev=[fb[0],*reversed(fb)];newdrev=[fb[0],*reversed(fp)]
 def cof(inv=False):
  for j in (reversed(range(width)) if inv else range(width)):
   bp,bb,beta,alpha=ff[j];ops=fake_undo_ops(bp,bb,beta,alpha)+common_padding_ops(bp,bb,[beta,alpha],False)
   if inv:controlled_linear(E.b,ops,flags[j],True)
   cofactor(m,E,bp,bb,beta,alpha,pinv,inv,True)
   if not inv:controlled_linear(E.b,ops,flags[j])
 if inverse:
  if tail:mobius(m,E,math,newrev,newdrev,pinv,tail,False,p_is_inverse=True)
  cof(True);schur(m,E,math,F,width,flags,True)
  mobius(m,E,math,oldrev,olddrev,pinv,F,True,p_is_inverse=True)
 else:
  mobius(m,E,math,oldrev,olddrev,pinv,F,False,p_is_inverse=True)
  schur(m,E,math,F,width,flags);cof()
  if tail:mobius(m,E,math,newrev,newdrev,pinv,tail,True,p_is_inverse=True)
 return fp,fb,C[:-2*width],pinv
