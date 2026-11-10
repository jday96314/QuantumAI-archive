"""A/D/p frame from two-slot lookahead heads and reversed quotients.

One terminal field word retains the otherwise lost normalization gauge.
This avoids symbolic resultants. Prefix correctness assumes every real
scalar pivot is nonzero; singular/completed inputs still invert literally.
No complete decoder or resource improvement is claimed by this module.
"""
from contextlib import contextmanager
import mixed_quotient_jacobi as convert
import generalized_gap_quotient as poly
import lookahead_jacobi as look
import shared_gap_quotient_peel as common
from reciprocal_inplace import apply as invert
from jacobi_decoder import cf

def zero_record(m,E,word,bit):
 with m.eq_test(E,word,0) as zero:E.cx(zero,bit)

def normalize(m,E,F,G,last_gap,inverse=False):
 """Normalize both leading coefficients, retaining D's original lead."""
 with m.opened(E,'quotient_frame_d_lead',12) as lead,m.opened(E,'quotient_frame_d_zero',1) as zero:
  def copy():return common.choose(m,E,G[-1],G[-2],last_gap,lead)
  def pad(back=False):
   indices=range(len(G)-2) if back else reversed(range(len(G)-2))
   for i in indices:E.swaps(last_gap,list(zip(G[i],G[i+1])))
  def top():E.swaps(last_gap,list(zip(G[-1],G[-2])))
  def unit():E.x(G[-1][0]);E.cx(zero[0],G[-1][0])
  def f_scale(back=False):
   with cf.armed(m,E,F[-1]),m.inv_word(E,F[-1]) as fi,m.opened(E,'quotient_frame_f_product',12) as tmp:
    for w in F[:-1]:m.mon_slot(E,w,F[-1] if back else fi,fi if back else F[-1],tmp)
  def g_scale(back=False):
   with cf.armed(m,E,lead),m.inv_word(E,lead) as di,m.opened(E,'quotient_frame_g_product',12) as tmp:
    if back:m.mon_slot(E,F[-1],di,lead,tmp)
    for w in G:m.mon_slot(E,w,lead if back else di,di if back else lead,tmp)
    if not back:m.mon_slot(E,F[-1],lead,di,tmp)
  if inverse:
   zero_record(m,E,G[-1],zero[0]);E.swaps(None,list(zip(lead,G[-1])));unit();pad(True);top();invert(m,E,F[-1]);g_scale(True);f_scale(True);zero_record(m,E,lead,zero[0]);E.free(copy())
  else:
   E.mac(copy());zero_record(m,E,lead,zero[0]);f_scale();g_scale();invert(m,E,F[-1]);top();pad();unit();E.swaps(None,list(zip(lead,G[-1])));zero_record(m,E,G[-1],zero[0])

@contextmanager
def parameters(m,E,W):
 n=len(W)//2;assert n%2==0 and n>=2
 with m.opened(E,'quotient_frame_gap_flags',n//2) as flags,m.opened(E,'quotient_frame_terminal',12) as c:
  E.x(c[0]);packets=[W[4*i:4*i+4] for i in range(n//2)]
  def preframe(back=False):
   items=list(zip(packets,flags))
   for H,flag in (reversed(items) if back else items):
    if not back:zero_record(m,E,H[0],flag)
    look.frame(m,E,H,flag,back)
    if back:zero_record(m,E,H[0],flag)
  def conversion(back=False):
   items=list(enumerate(zip(packets,flags)))
   for i,(H,flag) in (reversed(items) if back else items):convert.packet_convert(m,E,H,flag,W[4*i-2] if i else None,back)
  def canonical(back=False):
   items=list(zip(packets,flags))
   for H,flag in (reversed(items) if back else items):E.swaps(flag,list(zip(H[0],H[2])))
  def reverse(back=False):
   def all_pairs():
    for i in range(n//2):
     for a,b in zip(W[2*i:2*i+2],W[2*(n-1-i):2*(n-i)]):E.swaps(None,list(zip(a,b)))
   def correct(undo=False):
    items=list(enumerate(flags))
    for i,flag in (reversed(items) if undo else items):
     j=n-2-2*i
     for a,b in zip(W[2*j:2*j+2],W[2*j+2:2*j+4]):E.swaps(flag,list(zip(a,b)))
   if back:correct(True);all_pairs()
   else:all_pairs();correct()
  preframe();conversion();canonical();reverse()
  rev_gaps={n-2-2*i:flag for i,flag in enumerate(flags)};gaps=[rev_gaps.get(i) for i in range(n)];F,G=convert.base.layout(W,c)
  poly.synthesize(m,E,F,G,gaps);normalize(m,E,F,G,flags[-1]);yield [(),*reversed(F[:-1])],[(),*reversed(G[:-1])],F[-1]
  normalize(m,E,F,G,flags[-1],True);poly.synthesize(m,E,F,G,gaps,True);reverse(True);canonical(True);conversion(True);preframe(True);E.x(c[0])
