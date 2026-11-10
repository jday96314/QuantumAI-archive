"""Shared quotient reconstruction for normal and leading-gap two-slot packets.

The early prefix may contain ordinary Jacobi pairs. Later lookahead packets
encode (a,b,h,e), with a=0 selecting one degree-two quotient. Only the short
headers are selected; polynomial tail products are shared between charts.
An exact zero census excludes the imposed gap marker from completion degree.
Every operation has a literal inverse, including arbitrary singular inputs.
"""
import quotient_jacobi_synthesis as base
from reciprocal_inplace import apply as invert
from quantum_control import selected
import counter_quotient_jacobi as counter
import shared_gap_quotient_peel as shared
import lookahead_jacobi as look
cf=base.cf

def mon(m,E,word,scale,eraser,t,control=None):
 if control is None:m.mon_slot(E,word,scale,eraser,t)
 else:
  with selected(m,E,control) as (cm,ce):cm.mon_slot(ce,word,scale,eraser,t)

def pair_convert(m,E,a,b,prev=None,inverse=False,control=None):
 def multiply(back=False):
  if prev is not None:
   with cf.armed(m,E,prev),m.inv_word(E,prev) as pi,m.opened(E,'mixed_previous_product',12) as t:
    mon(m,E,a,pi if back else prev,prev if back else pi,t,control)
 def reciprocal():
  with cf.armed(m,E,a),m.opened(E,'mixed_head_product',12) as t:
   invert(m,E,a,control=control,before_swap=lambda src,iv:mon(m,E,b,iv,src,t,control))
 if inverse:reciprocal();multiply(True)
 else:multiply();reciprocal()

def packet_convert(m,E,H,flag,prev=None,inverse=False):
 a,b,c,d=H
 E.x(flag)
 if inverse:pair_convert(m,E,c,d,a,True,flag);pair_convert(m,E,a,b,prev,True,flag)
 else:pair_convert(m,E,a,b,prev,control=flag);pair_convert(m,E,c,d,a,control=flag)
 E.x(flag)
 def multiply(back=False):
  if prev is not None:
   with cf.armed(m,E,prev),m.inv_word(E,prev) as pi,m.opened(E,'mixed_gap_previous_product',12) as t:
    mon(m,E,b,pi if back else prev,prev if back else pi,t,flag)
 def reciprocal():
  with cf.armed(m,E,b),m.opened(E,'mixed_gap_head_product',12) as t:
   def scale(src,iv):
    for w in (c,d):mon(m,E,w,iv,src,t,flag)
   invert(m,E,b,control=flag,before_swap=scale)
 if inverse:E.swaps(flag,list(zip(b,c)));reciprocal();multiply(True)
 else:multiply();reciprocal();E.swaps(flag,list(zip(b,c)))

def inner_put(m,E,F,G,flag,inverse=False):
 """Normal inner linear fold; on gap leave both tail banks untouched."""
 q0,q1=G[-2:];D=F[:-2];R=G[:-2]
 def product():
  if not R:return
  E.x(flag)
  with shared.masks(m,E,[q0,q1],flag) as q:
   E.xor(m.polyprod_rt(q,D[:-1],[[r] for r in R]+[[]]))
  E.x(flag)
 def header(back=False):
  E.x(flag)
  with selected(m,E,flag) as (cm,ce),cf.armed(cm,ce,D[-1]),cm.inv_word(ce,D[-1]) as di,cm.opened(ce,'mixed_inner_product',12) as t:
   def scale(w):cm.mon_slot(ce,w,di if back else D[-1],D[-1] if back else di,t)
   if back:scale(q1)
   else:scale(q0)
   if len(D)>1:ce.xor(cm.prod_rt(q0,q1,D[-2]))
   if back:scale(q0)
   else:scale(q1)
  E.x(flag)
 if inverse:header(True);product()
 else:product();header()

def poly_put(m,E,F,G,packet_flags,inverse=False,index=0):
 if not G:return
 if len(G)==1:
  base.poly.peel(m,E,F,G,inverse=not inverse);return
 flag=packet_flags[index]
 def outer(back=False):
  if len(G)>=3:shared.peel(m,E,F,G,flag,back)
  else:
   # Common zero low coefficients extend the final two-slot packet to
   # the validated shared primitive's minimum nominal degree three.
   with m.opened(E,'mixed_final_low_f',12) as zf,m.opened(E,'mixed_final_low_g',12) as zg:
    shared.peel(m,E,[zf,*F],[zg,*G],flag,back)
 if inverse:
  if flag is None:
   base.poly.peel(m,E,F,G);base.poly.peel(m,E,G,F[:-2])
  else:
   outer()
   E.swaps(flag,list(zip(sum(F[:-2],[]),sum(G[:-1],[]))))
   inner_put(m,E,F,G,flag,True)
  poly_put(m,E,F[:-2],G[:-2],packet_flags,True,index+1)
 else:
  poly_put(m,E,F[:-2],G[:-2],packet_flags,False,index+1)
  if flag is None:
   base.poly.peel(m,E,G,F[:-2],True);base.poly.peel(m,E,F,G,True)
  else:
   inner_put(m,E,F,G,flag)
   E.swaps(flag,list(zip(sum(F[:-2],[]),sum(G[:-1],[]))))
   outer(True)

def synthesize(m,E,W,c,metadata,gap_flags,inverse=False,through=0):
 n=len(W)//2;count,armed=metadata[:-1],metadata[-1];F,G=base.layout(W,c)
 packets=[];gi=0
 for k in range(0,n-1,2):
  flag=None if k<through else gap_flags[gi]
  if flag is not None:gi+=1
  packets.append((k,W[2*k:2*k+4],flag))
 assert gi==len(gap_flags)
 def frame(back=False):
  for k,H,flag in (reversed(packets) if back else packets):
   if flag is not None:
    if not back:
     with m.eq_test(E,H[0],0) as z:E.cx(z,flag)
    look.frame(m,E,H,flag,back)
    if back:
     with m.eq_test(E,H[0],0) as z:E.cx(z,flag)
 def census(back=False):
  for k,H,flag in (reversed(packets) if back else packets):
   if flag is None:
    for a in (reversed(H[::2]) if back else H[::2]):
     with m.eq_test(E,a,0) as z:counter.increment(m,E,count,z,back)
   else:
    E.x(flag)
    with selected(m,E,flag) as (cm,ce):
     for a in (reversed(H[::2]) if back else H[::2]):
      with cm.eq_test(ce,a,0) as z:counter.increment(cm,ce,count,z,back)
    E.x(flag)
    with selected(m,E,flag) as (cm,ce),cm.eq_test(ce,H[1],0) as z:
     for _ in range(2):counter.increment(cm,ce,count,z,back)
  if n%2:
   with m.eq_test(E,W[-2],0) as z:counter.increment(m,E,count,z,back)
 def convert(back=False):
  if back and n%2:pair_convert(m,E,W[-2],W[-1],W[-4] if n>1 else None,True)
  for k,H,flag in (reversed(packets) if back else packets):
   prev=W[2*k-2] if k else None
   if flag is None:
    if back:pair_convert(m,E,H[2],H[3],H[0],True);pair_convert(m,E,H[0],H[1],prev,True)
    else:pair_convert(m,E,H[0],H[1],prev);pair_convert(m,E,H[2],H[3],H[0])
   else:packet_convert(m,E,H,flag,prev,back)
  if not back and n%2:pair_convert(m,E,W[-2],W[-1],W[-4] if n>1 else None)
 def route():
  for _,H,flag in packets:
   if flag is not None:E.swaps(flag,list(zip(H[0],H[2])))
 def terminal():
  with counter.equal(m,E,count,0) as full:E.cx(full,c[0])
 def record():
  for i,a in enumerate(W[::2]):
   with counter.equal(m,E,count,n-i) as match,m.eq_test(E,a,0) as z:E.b.ccx(match,z,armed)
 def flip():
  for i,a in enumerate(W[::2]):
   with counter.equal(m,E,count,n-i) as match:E.b.ccx(match,armed,a[0])
 flags=[f for _,_,f in packets]
 if inverse:
  poly_put(m,E,F,G,flags,True);flip();record();terminal();route();convert(True);census(True);frame(True)
 else:
  frame();census();convert();route();terminal();record();flip();poly_put(m,E,F,G,flags)
 return F,G
