"""Block head stored in denominator continuants plus one determinant scalar.

Nonzero Jacobi heads use D_j=(t+b_j)D_(j-1)+a_j D_(j-2). For a completed
zero suffix, coherently apply common-t padding instead. The head's 2b words
hold A's b coefficients, D's b-1 coefficients, and p; no extra A/D banks.
The b-1 zero flags are retained only while this reversible frame is active.
"""
from contextlib import contextmanager,ExitStack
from quantum_control import selected,permutation_ops,controlled_ops
from jacobi_decoder import cf
from field_dag import TreeEvaluator,Evaluator
from continuant_frame_affine import compiled


def layout(W):
 previous=[];current=[W[1]];frames=[]
 for i in range(1,len(W)//2):
  beta,alpha=W[2*i:2*i+2];frames.append((previous,current,beta,alpha));previous,current=current,previous+[beta,alpha]
 return frames,[(),*reversed(current)],[(),*reversed(previous)],W[0]


def frame(m,E,W,flags,inverse=False):
 frames,A,D,p=layout(W)
 def generic(cm,ce,P,B,beta,alpha,back=False):
  def high():
   ce.xor(cm.prod_rt(beta,alpha,B[-1]))
   if len(B)>1:
    for a,b in zip(B[-2],beta):ce.cx(a,b)
   for a,b in zip(B[-1],alpha):ce.cx(a,b)
  def undo_high():
   for a,b in zip(B[-1],alpha):ce.cx(a,b)
   if len(B)>1:
    for a,b in zip(B[-2],beta):ce.cx(a,b)
   ce.xor(cm.prod_rt(beta,alpha,B[-1]))
  if back:undo_high()
  with cm.inv_word(ce,beta) as bi,cm.opened(ce,'continuant_frame_product',12) as tmp:
   if not back:cm.mon_slot(ce,p,beta,bi,tmp)
   for i,w in enumerate(P):
    if not back:cm.mon_slot(ce,w,beta,bi,tmp)
    ce.xor(cm.prod_rt(w,alpha,B[i]))
    if i:
     for a,b in zip(B[i-1],w):ce.cx(a,b)
    if back:cm.mon_slot(ce,w,bi,beta,tmp)
   if back:cm.mon_slot(ce,p,bi,beta,tmp)
  if not back:high()
 for index in (reversed(range(len(frames))) if inverse else range(len(frames))):
  P,B,beta,alpha=frames[index];flag=flags[index];new=P+[beta,alpha]
  ops=permutation_ops([beta,*B,alpha,*P],[*new,*B])
  if not inverse:
   with m.eq_test(E,beta,0) as zero:E.cx(zero,flag)
  if inverse:controlled_ops(E,ops,flag,True)
  E.x(flag)
  with selected(m,E,flag) as (cm,ce):generic(cm,ce,P,B,beta,alpha,inverse)
  E.x(flag)
  if not inverse:controlled_ops(E,ops,flag)
  if inverse:
   with m.eq_test(E,beta,0) as zero:E.cx(zero,flag)
 return A,D,p


@contextmanager
def active_flag(m,E,flags):
 if len(flags)==1:
  E.x(flags[0]);yield flags[0];E.x(flags[0])
 elif not flags:
  with m.opened(E,'continuant_frame_active',1) as bit:
   E.x(bit[0]);yield bit[0];E.x(bit[0])
 else:
  with m.eq_test(E,flags,0) as active:yield active


def formula(m,E,math,c,root,variables,out,evaluator=TreeEvaluator,clear=False):
 ev=evaluator(m,E,math,c,variables)
 if not clear:ev.xor(root,out)
 else:
  with ev.value(root) as value:E.free([('L',tuple(sorted(s)),0,(t,)) for s,t in zip(value,out) if s])
 ev.clear()


class CachedEvaluator(Evaluator):
 def __init__(self,*args,**kwargs):super().__init__(*args,capacity=3,**kwargs)


@contextmanager
def parameters(m,E,math,H,evaluator=CachedEvaluator):
 b=len(H)//2
 with m.opened(E,'continuant_frame_zero_flags',b-1) as flags:
  A,D,p=frame(m,E,H,flags)
  with ExitStack() as stack:
   P=[stack.enter_context(m.opened(E,f'continuant_frame_affine_{j}',12)) for j in range(b)]
   @contextmanager
   def chart(d):
    with m.opened(E,'continuant_frame_chart',1) as word:
     bit=word[0]
     def toggle():
      if d==b:
       E.x(bit)
       if flags:E.cx(flags[-1],bit)
      else:
       E.cx(flags[d-1],bit)
       if d>1:E.cx(flags[d-2],bit)
     toggle();yield bit;toggle()
   def repay(rt,measurements):
    for term in rt:
     for t in term[-1]:
      cond=measurements[t]
      if term[0]=='L':
       for q in term[1]:E.b.z(q,cond=cond)
       if term[2]:E.b.neg(cond=cond)
      else:E.pivot_pair(term[1],term[2],term[3],term[4],lambda a,z,cond=cond:E.b.cz(a,z,cond=cond))
   def prepare(clear=False):
    measurements={q:E.b.hmr(q) for word in P for q in word} if clear else None
    for d in range(1,b+1):
     c,roots=compiled(d,b);variables=[*A[1:d+1],*D[1:d]]
     ev=evaluator(m,E,math,c,variables)
     with chart(d) as active,m.opened(E,'continuant_frame_resultant',12) as res:
      if d==1:E.x(res[0])
      else:ev.xor(roots[0],res)
      with cf.armed(m,E,res),m.opened(E,'continuant_frame_ratio',12) as ratio:
       def ratio_step(back=False):
        with m.inv_word(E,res) as ri,m.opened(E,'continuant_frame_active_scalar',12) as scalar:
         sr=[t for a,z in zip(p,scalar) for t in m.one_rt((active,),0,(a,),0,z)];E.mac(sr)
         rr=m.prod_rt(ratio,scalar,ri);(E.free if back else E.mac)(rr);E.free(sr)
       ratio_step()
       for root,out in zip(roots[1:],P):
        with ev.value(root) as val:
         rt=[];m.bil_terms(rt,m.KT,val,m.spec_plain(ratio),tgt=out)
         if clear:repay(rt,measurements)
         else:E.mac(rt)
       ratio_step(True)
      if d==1:E.x(res[0])
      else:
       with ev.value(roots[0]) as value:E.free([('L',tuple(sorted(v)),0,(z,)) for v,z in zip(value,res) if v])
     ev.clear()
   prepare();yield A,D,p,P;prepare(True)
  frame(m,E,H,flags,True)

# Computing the affine boundary before changing head coordinates avoids
# resultant divisions entirely. Keep the rational version for comparison.
parameters_resultant=parameters


@contextmanager
def parameters(m,E,math,H,evaluator=TreeEvaluator):
 from continuant_polynomials import derive_jacobi,Compiler
 from functools import lru_cache
 b=len(H)//2
 if not hasattr(parameters,'plans'):parameters.plans={}
 if b not in parameters.plans:
  c=Compiler(2*b);roots=[c.compile(p) for p in derive_jacobi(b)['P']];parameters.plans[b]=(c,roots)
 c,roots=parameters.plans[b]
 with ExitStack() as stack:
  P=[stack.enter_context(m.opened(E,f'continuant_direct_affine_{j}',12)) for j in range(b)]
  def affine(clear=False):
   ev=evaluator(m,E,math,c,H)
   measurements={q:E.b.hmr(q) for word in P for q in word} if clear else None
   for root,out in zip(roots,P):
    if not clear:ev.xor(root,out)
    else:
     with ev.value(root) as value:
      for source,target in zip(value,out):
       for q in source:E.b.z(q,cond=measurements[target])
   ev.clear()
  affine()
  with m.opened(E,'continuant_direct_zero_flags',b-1) as flags:
   A,D,p=frame(m,E,H,flags);yield A,D,p,P;frame(m,E,H,flags,True)
  affine(True)
