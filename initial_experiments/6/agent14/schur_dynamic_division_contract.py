"""Contract runtime masks into a single division phase field.

Create uses the compact in-place numerator from the preceding component.
Phase streams classical coefficient maps through the native bilinear MAC
instead of materializing both numerator fields and another contraction field.
"""
from functools import lru_cache
import kickmix as km
from schur_dynamic_division_compact import create,norm_rt
from extension_inverse import setup
from cubic_trace_field import inverse_columns,linear
from recursive_forney import inverse_dot_data
from tower_inverse import apply_map,map_ops


@lru_cache(None)
def maps():
 m,_,fw,bk,tr=setup();host=km.GF2Field(12)
 def nt(x):
  v=0
  for _ in range(12):v^=x;x=m.gmul(x,x)
  assert v in (0,1);return v
 def ht(x):
  v=0
  for _ in range(12):v^=x;x=host.mul(x,x)
  assert v in (0,1);return v
 tn=tuple(sum(nt(m.gmul(1<<i,1<<j))<<j for j in range(12)) for i in range(12));dual=inverse_columns(tn);th=tuple(sum(ht(host.mul(1<<i,1<<j))<<j for j in range(12)) for i in range(12));output=tuple(linear(th,tr(1<<i,bk)) for i in range(12))
 return dual,output


@lru_cache(None)
def columns(w):
 m=setup()[0];dual,out=maps();return tuple(linear(out,m.gmul(w,x)) for x in dual)


def phase(b,E,num,den,mask,kappa,controls,clear_completed=None):
 m,_,fw,bk,tr=setup();n0,n1=num[:12],num[12:];d0,d1=den[:12],den[12:]
 with m.opened(E,'contract_norm',12) as nn,m.opened(E,'contract_field',12) as uu,m.opened(E,'contract_bd',12) as bd:
  high=None if clear_completed is None else E.A.open('contract_high',12)
  nr=norm_rt(m,den,nn,kappa,controls);qr=m.prod_rt(bd,n1,d1);E.mac(nr);E.mac(qr);work=b.alloc_bits(2);co,parity=work
  def coefficient(w,words,action):
   cc=columns(w)
   for j in range(12):
    src=km.array([word[i] for word in words for i,c in enumerate(cc) if c>>j&1]);b.cx(src,co);action(j,co);b.cx(src,co)
  def classical_parity(w,words,measurement):
   b.bit_store(parity,False)
   def bit(j,control):b.push_condition(control);b.cx(measurement[j],parity);b.pop_condition()
   coefficient(w,words,bit)
  def pair(wa,ca,wb,cb,w,words,measurement=None,high_targets=()):
   if measurement is None:
    def mid(fan):
     def bit(j,control):E.b.cx(fan,uu[j],cond=control)
     coefficient(w,words,bit)
     for q in high_targets:E.b.cx(fan,q)
    E.fan_pair(wa,ca,wb,cb,mid)
   else:
    classical_parity(w,words,measurement);E.pivot_pair(wa,ca,wb,cb,lambda a,c:E.b.cz(a,c,cond=parity));b.bit_store(parity,False)
  def apply(measurement=None):
   rt=[]
   # Numerator low = P0 + kappa*P1, high = P0 + P2.
   m.bil_terms(rt,m.KT,m.spec_plain(n0),[a^c for a,c in zip(m.spec_plain(d0),m.spec_plain(d1))],tws_of=lambda js:tuple(js))
   for kind,wa,ca,wb,cb,targets in rt:
    pair(wa,ca,wb,cb,sum(1<<j for j in targets),[mask[:12],mask[12:]],measurement,() if high is None or measurement is not None else [high[j] for j in targets])
   rt=[];m.bil_terms(rt,m.KT,[a^c for a,c in zip(m.spec_plain(n0),m.spec_plain(n1))],m.spec_plain(d0),tws_of=lambda js:tuple(js))
   for kind,wa,ca,wb,cb,targets in rt:pair(wa,ca,wb,cb,sum(1<<j for j in targets),[mask[12:]],measurement,() if high is None or measurement is not None else [high[j] for j in targets])
   for i,q in enumerate(bd):
    w=m.gmul(kappa,1<<i)
    if measurement is None:coefficient(w,[mask[:12]],lambda j,control,q=q:E.b.cx(q,uu[j],cond=control))
    else:classical_parity(w,[mask[:12]],measurement);E.b.z(q,cond=parity);b.bit_store(parity,False)
    for flag,delta in controls:pair((flag,),0,(q,),0,m.gmul(delta,1<<i),[mask[:12]],measurement)
  apply()
  if high is not None:
   b.flip_if_equal(km.array([km.q(q) for q in high]),0,target=clear_completed)
   rt=[];m.bil_terms(rt,m.KT,m.spec_plain(n0),[a^c for a,c in zip(m.spec_plain(d0),m.spec_plain(d1))],high);m.bil_terms(rt,m.KT,[a^c for a,c in zip(m.spec_plain(n0),m.spec_plain(n1))],m.spec_plain(d0),high);E.free(rt);E.A.close('contract_high',high);high=None
  source=km.array([km.q(q) for q in nn]);target=km.array([km.q(q) for q in uu]);apply_map(b,source,map_ops(bk));leaf,inv,cleanup=inverse_dot_data(b,source,target);b.cz(leaf,inv);cleanup();apply_map(b,source,map_ops(bk),True)
  measurement=b.alloc_bits(12);b.hmr(target,measurement);apply(measurement);b.bit_store(measurement,False);b.free(measurement);b.bit_store(work,False);b.free(work);E.free(qr);E.free(nr)
