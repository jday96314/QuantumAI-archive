"""Dynamic quadratic extension inversion with coherent node controls."""
import kickmix as km
from extension_inverse import setup,init_inverse12,classical_transpose
from schur_dynamic_division import norm_rt
from tower_inverse import map_ops,apply_map
from recursive_forney import inverse_dot_data

def init_extension_inverse(m,E,source,out,kappa,controls):
 with m.opened(E,'extension_norm',12) as nn,m.opened(E,'extension_norm_inverse',12) as ni:
  nr=norm_rt(m,source,nn,kappa,controls);E.mac(nr)
  init_inverse12(m,E,nn,ni)
  left=[a^b for a,b in zip(m.spec_plain(source[:12]),m.spec_plain(source[12:]))]
  rt=[];m.bil_terms(rt,m.KT,left,m.spec_plain(ni),out[:12]);E.mac(rt)
  E.mac(m.prod_rt(out[12:],source[12:],ni))
  E.free(norm_rt(m,out,ni,kappa,controls));E.free(nr)

def inverse_phase(b,E,source,mask,kappa,controls):
 """mask dot inverse(source), with source stored in native tower coordinates."""
 m,_,fw,bk,tr=setup();a,c=source[:12],source[12:]
 with m.opened(E,'extension_phase_norm',12) as nn:
  rt=norm_rt(m,source,nn,kappa,controls);E.mac(rt)
  aa=km.array([km.q(i) for i in a]);cc=km.array([km.q(i) for i in c]);nr=km.array([km.q(i) for i in nn])
  # Native polynomial bits -> host polynomial bits. Classical masks use
  # the transpose of the opposite (host -> native) coordinate transform.
  fops=map_ops(fw);bops=map_ops(bk);trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans)
  for word in [aa,cc,nr]:apply_map(b,word,bops)
  for word in [mask[:12],mask[12:]]:apply_map(b,word,mops)
  u=b.alloc_qubits(12)
  def umap():
   b.cx(cc,aa);classical_transpose(b,aa,mask[:12],u);b.cx(cc,aa)
   classical_transpose(b,cc,mask[12:],u)
  umap();leaf,inv,cleanup=inverse_dot_data(b,nr,u);b.cz(leaf,inv);cleanup();umap();b.free(u)
  for word in [mask[:12],mask[12:]]:apply_map(b,word,mops,True)
  for word in [aa,cc,nr]:apply_map(b,word,bops,True)
  E.free(rt)

