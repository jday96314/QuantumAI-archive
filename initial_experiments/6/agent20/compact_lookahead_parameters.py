"""In-head Mobius parameters for ordinary/leading-zero two-pair blocks.

Encoded H=(a,b,h,e). Prepare the boundary before reframing. Put
u=b*(a+flag), v=a*h, k=h+e^2. Then P0=v*k+u*e^3 and
P1=h*(v*e+u*k), valid on both charts. The head itself stores A=(1,e,h),
D=(1,b) and p after a*=armed(h+b*e+b^2), followed by a/b exchange on the gap.
A zero ordinary c uses an invertible scalar1; completed residuals vanish,
while unsupported second-pivot gaps retain a defined literal inverse.
"""
from contextlib import contextmanager
from jacobi_decoder import cf


def repay(E,terms,measurements):
 for t in terms:
  for target in t[-1]:
   cond=measurements[target]
   if t[0]=='L':
    for q in t[1]:E.b.z(q,cond=cond)
    if t[2]:E.b.neg(cond=cond)
   else:E.pivot_pair(t[1],t[2],t[3],t[4],lambda a,b,cond=cond:E.b.cz(a,b,cond=cond))


def reframe(m,E,math,H,flag,inverse=False):
 a,b,h,e=H
 if inverse:E.swaps(flag,list(zip(a,b)))
 with m.opened(E,'compact_lookahead_c',12) as c:
  cr=m.prod_rt(c,b,e)+[('L',(q,),0,(z,)) for q,z in zip(h,c)]
  m.lin_rows_terms(cr,math.rows(1,2),b,c);E.mac(cr)
  with cf.armed(m,E,c),m.inv_word(E,c) as ci,m.opened(E,'compact_lookahead_scale',12) as temp:
   m.mon_slot(E,a,ci if inverse else c,c if inverse else ci,temp)
  E.free(cr)
 if not inverse:E.swaps(flag,list(zip(a,b)))


@contextmanager
def parameters(m,E,math,H,flag):
 a,b,h,e=H
 with m.opened(E,'compact_lookahead_P0',12) as p0,m.opened(E,'compact_lookahead_P1',12) as p1:
  def boundary(clear=False):
   measurements={q:E.b.hmr(q) for w in (p0,p1) for q in w} if clear else None
   with m.opened(E,'compact_lookahead_u',12) as u,m.opened(E,'compact_lookahead_v',12) as v:
    ap=m.spec_plain(a);ap[0]=frozenset(ap[0])^frozenset((flag,));ur=[];m.bil_terms(ur,m.KT,m.spec_plain(b),ap,tgt=u)
    vr=m.prod_rt(v,a,h);E.mac(ur);E.mac(vr)
    k=[x^y for x,y in zip(m.spec_plain(h),math.spec([(1,e,2)]))]
    with m.opened(E,'compact_lookahead_cube',12) as cube:
     er=math.cube_rts(m.spec_plain(e),lambda bits:[cube[j] for j in bits]);E.mac(er)
     rt=m.prod_rt(p0,u,cube);m.bil_terms(rt,m.KT,m.spec_plain(v),k,tgt=p0)
     if clear:repay(E,rt,measurements)
     else:E.mac(rt)
     E.free(er)
    with m.opened(E,'compact_lookahead_P1_inner',12) as inner:
     ir=m.prod_rt(inner,v,e);m.bil_terms(ir,m.KT,m.spec_plain(u),k,tgt=inner);E.mac(ir)
     rt=m.prod_rt(p1,h,inner)
     if clear:repay(E,rt,measurements)
     else:E.mac(rt)
     E.free(ir)
    E.free(vr);E.free(ur)
  boundary();reframe(m,E,math,H,flag);yield [(),e,h],[(),b],a,(p0,p1);reframe(m,E,math,H,flag,True);boundary(True)
