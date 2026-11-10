"""Gap-specialized three-degree compact block with exact tail reconstruction.

Domain: the local series has first moment zero and next pivot and degree-three
normalizer nonzero. The first moment is retained as an exceptional coordinate.
This is a reversible branch component, without global branch dispatch.
"""
from gap_cofactor import apply as cof
from mobius_series import apply as mobius
from tower_inverse import map_ops

def head(m,E,F,inverse=False):
 e,d,a,b,c,f=F
 def square_xor(src,dst):
  for i,q in enumerate(src):
   v=m.gmul(1<<i,1<<i)
   for j in range(12):
    if v>>j&1:E.cx(q,dst[j])
 def square(word,back=False):
  ops=map_ops([m.gmul(1<<i,1<<i) for i in range(12)])
  for kind,i,j in (reversed(ops) if back else ops):getattr(E.b,kind)(word[i],word[j])
 def cube():
  square(a);E.xor(m.prod_rt(c,a,F[0])) if False else None
  # a is now squared; recover unsquared a in a clean temporary for the cube.
  with m.opened(E,'gap_unsquared_a',12) as tmp:
   for q,t in zip(a,tmp):E.cx(q,t)
   square(tmp,True);E.xor(m.prod_rt(c,a,tmp))
   square(tmp)
   for q,t in zip(a,tmp):E.cx(q,t)
  square(a,True)
 def last_shear():
  square_xor(b,f);square(a);E.xor(m.prod_rt(f,b,a));square_xor(a,f);square(a,True)
 def scale_d(back):
  with m.inv_word(E,d) as di,m.opened(E,'gap_head_scale',12) as t:
   for word in [a,b,c,f]:m.mon_slot(E,word,d if back else di,di if back else d,t)
 def scale_gamma(back):
  with m.inv_word(E,c) as ci,m.opened(E,'gap_head_scale',12) as t:m.mon_slot(E,f,c if back else ci,ci if back else c,t)
 if inverse:scale_gamma(True);last_shear();cube();square_xor(a,b);scale_d(True)
 else:scale_d(False);square_xor(a,b);cube();last_shear();scale_gamma(False)
 return [d,b,a,c,f,e]

def apply(m,E,math,Bp,B,C,pinv,inverse=False):
 physical=list(reversed(C))[:6];tail=list(reversed(C))[6:];rev=[B[0],*reversed(B)];drev=[B[0],*reversed(Bp)]
 def local(T,back=False):
  if back:head(m,E,T,True)
  mobius(m,E,math,rev,drev,pinv,T,back,p_is_inverse=True)
  if not back:head(m,E,T)
 with m.opened(E,'gap_block_head',72) as bank:
  raw=[bank[12*i:12*(i+1)] for i in range(6)]
  q=[raw[i] for i in [1,3,2,4,5,0]]
  fp=Bp+q[:3];fb=B+[q[3],q[5],q[4]];finalrev=[fb[0],*reversed(fb)]
  def copy():
   for source,target in zip(physical,raw):
    for x,y in zip(source,target):E.cx(x,y)
  def swap():
   for a,b in zip(raw,physical):E.swaps(None,list(zip(a,b)))
  if not inverse:
   copy();local(raw);cof(m,E,Bp,B,q,pinv,norm=False)
   math.triangular_suffix(finalrev,1,list(reversed(C)),6)
   cof(m,E,Bp,B,q,pinv,True,norm=False)
   math.triangular(rev,1,tail,True)
   local(physical)
   for source,target in zip(raw,physical):
    for x,y in zip(source,target):E.cx(x,y)
   cof(m,E,Bp,B,q,pinv,norm=True);swap()
  else:
   swap();cof(m,E,Bp,B,q,pinv,True,norm=True)
   for source,target in zip(raw,physical):
    for x,y in zip(source,target):E.cx(x,y)
   local(physical,True);math.triangular(rev,1,tail)
   cof(m,E,Bp,B,q,pinv,norm=False)
   math.triangular_suffix(finalrev,1,list(reversed(C)),6,True)
   cof(m,E,Bp,B,q,pinv,True,norm=False);local(raw,True);copy()
  sub={a:b for src,dst in zip(raw,physical) for a,b in zip(src,dst)}
  convert=lambda words:[[sub.get(q,q) for q in word] for word in words]
  return convert(fp),convert(fb),C[:-6],pinv
