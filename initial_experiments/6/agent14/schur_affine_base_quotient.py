"""Recover the two heads from a degree-one affine Schur output."""
import kickmix as km
from extension_inverse import setup
from low_workspace_inverse import bind


def create(b,E,n0,n1,b0,u,t,kappa,controls):
 m=setup()[0];m.inv_word=bind(m);U=[q.id for q in u];T=[q.id for q in t];N0=[q.id for q in n0];N1=[q.id for q in n1];B0=[q.id for q in b0]
 def hmap():
  for i,q in enumerate(b0):
   c=m.gmul(1<<i,1<<i)^(1<<i)
   for j in range(12):
    if c>>j&1:b.cx(q,u[j])
  for j in range(12):
   if kappa>>j&1:b.x(u[j])
   for flag,delta in controls:
    if delta>>j&1:b.cx(km.q(flag),u[j])
 hmap();b.cx(n0,t);E.mac(m.prod_rt(T,N1,B0))
 with m.inv_word(E,U) as hi,m.opened(E,'base_quotient',12) as tmp:m.mon_slot(E,T,hi,U,tmp)
 hmap();b.cx(n1,u);b.cx(t,u);E.mac(m.prod_rt(U,T,B0))
