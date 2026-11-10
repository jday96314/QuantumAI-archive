"""Specialize products x*x^2 to a verified 28-AND cubing cover."""
import json
from pathlib import Path
from functools import lru_cache
from triangular_dp import DPArithmetic
from series_inverse import SeriesArithmetic

@lru_cache(None)
def cube_data():return json.loads(Path(__file__).with_name('parameters').joinpath('cube28.json').read_text())

class GoldArithmetic(DPArithmetic):
 def frobenius(self,v,q):
  rows=self.rows(1,q);out=[]
  for mask in rows:
   s=set()
   for j in self.m.bitlist(mask,12):s.symmetric_difference_update(v[j])
   out.append(frozenset(s))
  return out
 def cube_terms(self,v,targets):
  z=cube_data();rt=[];m=self.m
  def out(mask):return targets(m.bitlist(mask,12))
  for a,b,c in z['terms']:
   term=m.norm_term(m.fold(v,m.bitlist(a,12)),0,m.fold(v,m.bitlist(b,12)),0,out(c))
   if term is not None:rt.append(term)
  for i,c in enumerate(z['linear_columns']):
   dst=out(c)
   if c and v[i] and dst:rt.append(('L',tuple(sorted(v[i])),0,dst))
  return rt
 def cube(self,v,targets):
  self.E.emit_mac(self.cube_terms(v,targets))
 def product(self,left,right,targets):
  source=None
  if list(right)==self.frobenius(left,2):source=left
  elif list(left)==self.frobenius(right,2):source=right
  if source is None:return SeriesArithmetic.product(self,left,right,targets)
  generic=[]
  self.m.bil_terms(generic,self.m.KT,left,right,tws_of=targets)
  if hasattr(self.m,'branch_flag'):
   base=self.m.base_module;E=self.m.base_emitter
   if hasattr(base,'branch_flag'):return self.E.emit_mac(generic)
   pp=sum(t[0]=='P' for t in generic);ll=len(generic)-pp
   right_bits=set(q for t in generic if t[0]=='P' for q in t[3])
   cost=(pp+len(right_bits)) if not ll and len(right_bits)<=12 else 2*pp+ll
   if cost<=40 or E.b.active_qubits+12>1004:return self.E.emit_mac(generic)
   with base.opened(E,'controlled_gold_source',12) as masked:
    rt=[t for i,w in enumerate(masked) for t in base.one_rt((self.m.branch_flag,),0,tuple(sorted(source[i])),0,w)]
    E.mac(rt)
    GoldArithmetic(base,E).cube(GoldArithmetic(base,E).spec([(1,masked,1)]),targets)
    E.free(rt)
   return
  specialized=self.cube_terms(source,targets)
  price=lambda rt:sum(t[0]=='P' for t in rt)
  self.E.emit_mac(specialized if price(specialized)<price(generic) else generic)
