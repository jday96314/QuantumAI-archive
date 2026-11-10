"""Fixed-slot polynomial continuants with a linear, quadratic or skipped slot.

A quadratic quotient occupies its slot and the following skipped slot.
Both real charts share one long remainder product; skipped slots mask all
coefficients to zero and exchange the tail banks. No chart controls a whole
long polynomial multiplication. Overlapping arbitrary flags give skip
priority and remain a literal reversible map.
"""
from contextlib import contextmanager
import quotient_jacobi_synthesis as base
from jacobi_decoder import cf
import shared_gap_quotient_peel as common

@contextmanager
def charts(m,E,gap,skip):
 with m.opened(E,'general_quotient_charts',2) as bank:
  effective,normal=bank
  def prepare(back=False):
   def ge():
    if gap is None:return
    if skip is None:E.cx(gap,effective)
    else:
     E.x(skip);E.b.ccx(skip,gap,effective);E.x(skip)
   def no():
    E.x(normal);E.cx(effective,normal)
    if skip is not None:E.cx(skip,normal)
   if back:no();ge()
   else:ge();no()
  prepare();yield effective,normal;prepare(True)

def choice(m,E,normal,gap,gap_flag,skip,target,one_off=False):
 rt=common.linear(m,E,[normal],target)
 for i,(a,b,t) in enumerate(zip(normal,gap,target)):
  rt+=m.one_rt((gap_flag,),0,(a,b),0,t)
  if skip is not None:rt+=m.one_rt((skip,),0,(a,),int(one_off and i==0),t)
 return rt

def peel(m,E,F,G,gap,skip=None,inverse=False):
 n=len(G);assert n>=3 and len(F)==n+1
 with charts(m,E,gap,skip) as (g,normal):
  def quotient(back=False):
   with m.opened(E,'general_quotient_lead',12) as lead:
    lr=choice(m,E,G[-1],G[-2],g,skip,lead,True);E.mac(lr)
    with cf.armed(m,E,lead),m.inv_word(E,lead) as iv,m.opened(E,'general_quotient_product',12) as tmp:
     def scale(word,reverse=False):m.mon_slot(E,word,lead if reverse else iv,iv if reverse else lead,tmp)
     def third_shear():
      with m.opened(E,'general_quotient_shear_sum',12) as out:
       rt=m.prod_rt(out,F[-2],G[-3])
       if n>=4:rt+=m.prod_rt(out,F[-1],G[-4])
       E.mac(rt)
       for a,b in zip(out,F[-3]):E.b.ccx(g,a,b)
       E.free(rt)
     def second_shear():
      with m.opened(E,'general_quotient_sublead',12) as sub:
       sr=choice(m,E,G[-2],G[-3],g,skip,sub);E.mac(sr);E.xor(m.prod_rt(F[-2],F[-1],sub));E.free(sr)
     def third_scale(reverse=False):
      with m.opened(E,'general_quotient_routed',12) as routed:
       E.swaps(g,list(zip(F[-3],routed)));m.mon_slot(E,routed,lead if reverse else iv,iv if reverse else lead,tmp);E.swaps(g,list(zip(F[-3],routed)))
     if back:third_scale(True);third_shear();scale(F[-2],True);second_shear();scale(F[-1],True)
     else:scale(F[-1]);second_shear();scale(F[-2]);third_shear();third_scale()
    E.free(lr)
  def remainder():
   with m.opened(E,'general_quotient_coefficients',36) as bank,m.opened(E,'general_quotient_boundary',12) as boundary:
    q0,q1,q2=[bank[12*i:12*(i+1)] for i in range(3)];rt=choice(m,E,F[-2],F[-3],g,skip,q0)+choice(m,E,F[-1],F[-2],g,skip,q1)
    for a,b in zip(F[-1],q2):rt+=m.one_rt((g,),0,(a,),0,b)
    E.mac(rt);cr=m.polyprod_rt([q0,q1,q2],G[:-1],[[r] for r in F[:-3]]+[[boundary],[],[]]);E.xor(cr)
    for a,b in zip(boundary,F[-3]):E.b.ccx(normal,a,b)
    bs=set(boundary);clear=[(*t[:-1],tuple(q for q in t[-1] if q in bs)) for t in cr];E.free([t for t in clear if t[-1]]);E.free(rt)
  if inverse:remainder();quotient(True)
  else:quotient();remainder()

def final(m,E,F,G,skip=None,inverse=False):
 assert len(F)==2 and len(G)==1
 with m.opened(E,'general_final_lead',12) as lead:
  rt=common.linear(m,E,G,lead)
  if skip is not None:
   for i,(a,b) in enumerate(zip(G[0],lead)):rt+=m.one_rt((skip,),0,(a,),int(i==0),b)
  E.mac(rt)
  with cf.armed(m,E,lead),m.inv_word(E,lead) as iv,m.opened(E,'general_final_product',12) as tmp:
   for word in (F if inverse else reversed(F)):m.mon_slot(E,word,lead if inverse else iv,iv if inverse else lead,tmp)
  E.free(rt)

def last_two(m,E,F,G,gap,skip=None,inverse=False):
 """Final linear/quadratic/skip step without padded zero field words."""
 assert len(F)==3 and len(G)==2
 with charts(m,E,gap,skip) as (g,normal),m.opened(E,'general_last_two_lead',12) as lead:
  lr=choice(m,E,G[-1],G[-2],g,skip,lead,True);E.mac(lr)
  def quotient(back=False):
   with cf.armed(m,E,lead),m.inv_word(E,lead) as iv,m.opened(E,'general_last_two_product',12) as tmp:
    def scale(word):m.mon_slot(E,word,lead if back else iv,iv if back else lead,tmp)
    def shear():
     with common.masks(m,E,[F[-1]],normal) as q:E.xor(m.prod_rt(F[-2],q[0],G[0]))
    def third():
     with m.opened(E,'general_last_two_routed',12) as routed:
      E.swaps(g,list(zip(F[0],routed)));m.mon_slot(E,routed,lead if back else iv,iv if back else lead,tmp);E.swaps(g,list(zip(F[0],routed)))
    if back:third();scale(F[-2]);shear();scale(F[-1])
    else:scale(F[-1]);shear();scale(F[-2]);third()
  def remainder():
   with common.masks(m,E,[F[-2]],normal) as q:E.xor(m.prod_rt(F[0],q[0],G[0]))
  if inverse:remainder();quotient(True)
  else:quotient();remainder()
  E.free(lr)

def synthesize(m,E,F,G,gaps,inverse=False,index=0):
 if not G:return
 skip=gaps[index-1] if index else None;gap=gaps[index]
 def step(back=False):
  if len(G)==1:final(m,E,F,G,skip,back)
  elif len(G)==2:last_two(m,E,F,G,gap,skip,back)
  else:peel(m,E,F,G,gap,skip,back)
 def route():
  if skip is not None:E.swaps(skip,list(zip(sum(G,[]),sum(F[:-1],[]))))
 if inverse:step();route();synthesize(m,E,G,F[:-2],gaps,True,index+1)
 else:synthesize(m,E,G,F[:-2],gaps,False,index+1);route();step(True)
