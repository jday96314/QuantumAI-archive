"""Outer S0 D0 fold when a0/b0 is the old numerator leading field.

The endpoint ratio is carried in the original a0 slot through the inner
chain. This fold consumes only b0 and drops its zero numerator leading slot.
It retains the same rare identity/completion chart convention as canonical.
"""
import kickmix as km
from extension_inverse import Adapter,setup,classical_transpose
from schur_dynamic_division_compact import numerator,norm_rt
from schur_dynamic_division_contract import phase as division_phase
from schur_dynamic_fold_recycle import evaluate,dapply
from low_workspace_inverse import bind
from schur_control import caxpy
from tower_inverse import apply_map,map_ops
from recursive_forney import inverse_dot_data
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create
from schur_dynamic_division_sequential import create as sequential_create


def ratio_phase(b,E,low,high,mask):
 m,_,fw,bk,tr=setup();a=km.array([km.q(q) for q in low]);v=km.array([km.q(q) for q in high]);ops=map_ops(bk);trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans)
 apply_map(b,a,ops);apply_map(b,v,ops);apply_map(b,mask,mops);u=b.alloc_qubits(12);classical_transpose(b,a,mask,u);leaf,iv,cleanup=inverse_dot_data(b,v,u);b.cz(leaf,iv);cleanup();classical_transpose(b,a,mask,u);b.free(u);apply_map(b,mask,mops,True);apply_map(b,v,ops,True);apply_map(b,a,ops,True)


def build(d,kappa,deltas,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());headers=[];rid=0
 def reg(n,name):
  nonlocal rid
  v=b.alloc_qubits(n)
  for q in v:headers.append(f'APPEND_TO_REGISTER q{q.id} r{rid}')
  headers.append(f'REGISTER r{rid} "{name}"');rid+=1;return v
 bv=reg(12,'b');N=[reg(12,f'N{i}') for i in range(d+1)];B=[reg(12,f'B{i}') for i in range(d)];nl=reg(12,'Nlead');bt=reg(12,'Btop');flags=reg(len(deltas),'flags');controls=[(q.id,v) for q,v in zip(flags,deltas)];C=b.alloc_qubits(1)[0];ghost=b.alloc_bits(25);headers.extend(f'APPEND_TO_REGISTER b{v.id} r{rid}' for v in ghost);headers.append(f'REGISTER r{rid} "head_ghost"');ab=Adapter(b);sc,E=lean_create(m,ab);BB=B+[bt]
 def classify():b.flip_if_equal(km.array(list(N[-1])+list(bv)),1<<12,target=C)
 def correction():
  for bit in range(12):
   if kappa>>bit&1:b.cx(C,N[-1][bit])
  b.cx(C,bt[0]);caxpy(m,E,C.id,[q.id for q in nl],[q.id for q in bt],1)
 def affine(back=False):
  def scale():
   T=[q.id for q in bv]
   with m.inv_word(E,T) as ti,m.opened(E,'outer_scale',12) as tmp:
    for n in N:m.mon_slot(E,[q.id for q in n],ti if back else T,T if back else ti,tmp)
  def shear():
   for n,v in zip(N,BB):E.mac(m.prod_rt([q.id for q in n],[q.id for q in nl],[q.id for q in v]))
  if back:scale();shear()
  else:shear();scale()
 if inverse:
  with evaluate(b,E,N,BB,kappa,controls) as (nv,dv):
   out=[q.id for q in list(nl)+list(bv)]
   sequential_create(b,E,nv,dv,out,kappa,controls);b.flip_if_equal(bv,0,target=C);b.cx(C,bv[0]);T=[q.id for q in bv]
   with m.inv_word(E,T) as ti,m.opened(E,'outer_ratio',12) as tmp:m.mon_slot(E,[q.id for q in nl],ti,T,tmp)
  for q,bit in zip(list(nl)+list(bv),ghost[:24]):ab.z(q.id,bit)
  ab.z(C.id,ghost[24]);b.bit_store(ghost,False)
  affine(True);correction();dapply(b,E,N,BB,kappa,controls,True);b.x(bt[0]);b.cx(N[-1],nl);classify()
 else:
  b.bit_store(ghost,False);classify();b.cx(N[-1],nl);b.x(bt[0]);dapply(b,E,N,BB,kappa,controls);correction();affine();b.hmr(km.array(list(nl)+list(bv)),ghost[:24]);b.free(nl);b.free(bv);b.hmr(C,ghost[24]);b.free(km.array([C]))
 if inverse:
  reset=b.alloc_bits(1)[0];b.hmr(C,reset);b.bit_store(reset,False);b.free(km.array([reset]));b.free(km.array([C]))
 ab.release_recycle('outer_emitter',sc);c=b.finish_circuit();return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(25,c.num_bits)))
