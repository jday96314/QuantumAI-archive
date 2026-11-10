"""One shared polynomial tail product for degree-one and degree-two quotients.

F has nominal degree n, G degree n-1. On the gap chart G[-1]=0 and
G[-2]!=0, so the quotient has degree two. The quotient is stored in the
last two/three F words, respectively. One selected three-word quotient and
one common tail convolution replace duplicate chart-controlled long products.
The input flag stays live; both directions are literal inverses for arbitrary
flags and armed zero pivots. This primitive does not schedule whole Euclid.
"""
from contextlib import contextmanager
from jacobi_decoder import cf


def linear(m,E,words,target):
 return [('L',(q,),0,(t,)) for word in words for q,t in zip(word,target)]


def choose(m,E,normal,gap,flag,target):
 rt=linear(m,E,[normal],target)
 for a,b,t in zip(normal,gap,target):rt+=m.one_rt((flag,),0,(a,b),0,t)
 return rt


@contextmanager
def masks(m,E,words,flag):
 with m.opened(E,'gap_quotient_masks',12*len(words)) as bank:
  out=[bank[12*i:12*(i+1)] for i in range(len(words))];rt=[]
  for src,dst in zip(words,out):
   for a,b in zip(src,dst):rt+=m.one_rt((flag,),0,(a,),0,b)
  E.mac(rt);yield out;E.free(rt)


def peel(m,E,F,G,flag,inverse=False):
 n=len(G)
 if len(F)!=n+1 or n<3:raise ValueError('requires degree n>=3 and nominal degree n-1 G')
 def quotient(back=False):
  with m.opened(E,'shared_gap_lead',12) as lead,m.opened(E,'shared_gap_sublead',12) as sub:
   lr=choose(m,E,G[-1],G[-2],flag,lead);sr=choose(m,E,G[-2],G[-3],flag,sub)
   E.mac(lr);E.mac(sr)
   with cf.armed(m,E,lead),m.inv_word(E,lead) as iv,m.opened(E,'shared_gap_scale_product',12) as temp:
    def scale(word,rev=False):m.mon_slot(E,word,lead if rev else iv,iv if rev else lead,temp)
    def third_shear():
     with masks(m,E,[F[-1],F[-2]],flag) as (hi,mid):
      rt=m.prod_rt(F[-3],mid,G[-3])
      if n>=4:rt+=m.prod_rt(F[-3],hi,G[-4])
      E.xor(rt)
    def third_scale(rev=False):
     with m.opened(E,'shared_gap_routed_scale',12) as routed:
      E.swaps(flag,list(zip(F[-3],routed)))
      m.mon_slot(E,routed,lead if rev else iv,iv if rev else lead,temp)
      E.swaps(flag,list(zip(F[-3],routed)))
    if back:
     third_scale(True);third_shear();scale(F[-2],True);E.xor(m.prod_rt(F[-2],F[-1],sub));scale(F[-1],True)
    else:
     scale(F[-1]);E.xor(m.prod_rt(F[-2],F[-1],sub));scale(F[-2]);third_shear();third_scale()
   E.free(sr);E.free(lr)
 def remainder():
  with m.opened(E,'shared_gap_q_coefficients',36) as bank,m.opened(E,'shared_gap_remainder_boundary',12) as boundary:
   q0,q1,q2=[bank[12*i:12*(i+1)] for i in range(3)]
   rt=choose(m,E,F[-2],F[-3],flag,q0)+choose(m,E,F[-1],F[-2],flag,q1)
   for a,b in zip(F[-1],q2):rt+=m.one_rt((flag,),0,(a,),0,b)
   E.mac(rt)
   # F[:-3] is the common low remainder. The normal chart has one
   # additional coefficient, accumulated from a temporary field word.
   targets=[[r] for r in F[:-3]]+[[boundary],[],[]]
   cr=m.polyprod_rt([q0,q1,q2],G[:-1],targets);E.xor(cr)
   E.x(flag)
   for a,b in zip(boundary,F[-3]):E.b.ccx(flag,a,b)
   E.x(flag)
   # Boundary is still a function of stable selected q and G. Erase only
   # its targets; the low remainder is the actual output, not garbage.
   bs=set(boundary);clear=[(*t[:-1],tuple(x for x in t[-1] if x in bs)) for t in cr]
   E.free([t for t in clear if t[-1]])
   E.free(rt)
 if inverse:remainder();quotient(True)
 else:quotient();remainder()
 return G,F[:-2],(F[-2],F[-1])
