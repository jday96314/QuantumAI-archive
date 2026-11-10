"""Reversible fixed-node Schur frames through isolated scalar pivots.

Zero b remains in the field bank as a chart tag. Forward and inverse
use the previous two field pairs instead of a separate node permutation.
"""
from contextlib import contextmanager
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from base_inverse import inplace as invert12
from tower_inverse import map_ops,apply_map
from schur_control import caxpy,caxpy_rt
from schur_fold import axpy
from fast_dqi import FastBuilder


def build(node,inverse=False,boundary=2):
 m,k,fw,bk,tr=setup();m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());lo=b.create_quantum_register(12,'a');hi=b.create_quantum_register(12,'b');g=b.create_quantum_register(12,'g');pa=b.create_quantum_register(12,'prev_a');pb=b.create_quantum_register(12,'prev_b');ppa=b.create_quantum_register(12,'prior_a');ppb=b.create_quantum_register(12,'prior_b');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));A=[q.id for q in lo];B=[q.id for q in hi];G=[q.id for q in g];PB=[q.id for q in pb]
 @contextmanager
 def prevzero():
  if boundary==0:yield None
  else:
   with m.eq_test(E,PB,0) as z:yield z
 def offset(zp):
  if boundary==0:return
  b.cx(pa,lo);b.cx(pb,lo);b.cx(pb,hi);caxpy(m,E,zp,[q.id for q in pa],A,1)
  if boundary>1:
   with m.eq_test(E,[q.id for q in ppb],0) as zz,m.opened(E,'gap_previous_chart',1) as flag:
    E.x(zp);rr=[('P',(zp,),0,(zz,),0,(flag[0],))];E.mac(rr);E.x(zp)
    b.cx(pb,pa);b.cx(ppa,pa);caxpy(m,E,flag[0],[q.id for q in pa],A,1);b.cx(ppa,pa);b.cx(pb,pa);caxpy(m,E,flag[0],PB,B,1)
    E.x(zp);E.free(rr);E.x(zp)
 def scale(back=False):
  with m.inv_word(E,G) as gi,m.opened(E,'gap_scale_values',12) as temp:
   for word in [A,B]:m.mon_slot(E,word,gi if back else G,G if back else gi,temp)
 def gupdate(zp):
  invert12(b,g)
  with m.eq_test(E,B,0) as zc:
   E.cx(zc,B[0]);ops=map_ops([m.gmul(1<<i,1<<i) for i in range(12)]);apply_map(b,hi,ops)
   if zp is None:
    with m.inv_word(E,B) as bi,m.opened(E,'gap_scale_frame',12) as temp:m.mon_slot(E,G,B,bi,temp)
   else:
    with m.opened(E,'gap_frame_factor',12) as factor:
     rt=[('L',(q,),0,(t,)) for q,t in zip(B,factor)]+caxpy_rt(m,zp,B,factor,1)+[('L',(zp,),0,(factor[0],))];E.mac(rt)
     with m.inv_word(E,factor) as fi,m.opened(E,'gap_scale_frame',12) as temp:m.mon_slot(E,G,factor,fi,temp)
     E.free(rt)
   apply_map(b,hi,ops,True);E.cx(zc,B[0])
 with prevzero() as zp:
  if inverse:gupdate(zp);offset(zp);scale(True);axpy(b,hi,lo,node,m)
  else:axpy(b,hi,lo,node,m);scale();offset(zp);gupdate(zp)
 ab.release_recycle('gap_pivot_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(shots=2048):
 m,k,fw,bk,tr=setup();node=tr(26,fw);rng=np.random.default_rng(142059);rows=[]
 for boundary in [0,1,2]:
  aa,bb,pa,pb,ppa,ppb=rng.integers(0,4096,(6,shots));gg=rng.integers(1,4096,shots);pb[::5]=0;ppb[::7]=0
  if boundary==0:A=np.zeros(shots,dtype=int);B=A.copy()
  else:
   A=np.where(pb==0,0,pa^pb);B=pb.copy()
   if boundary>1:A=np.where((pb!=0)&(ppb==0),ppa,A);B=np.where((pb!=0)&(ppb==0),0,B)
  bb[::11]=np.array([m.gmul(int(p),m.ginv(int(g))) for p,g in zip(B[::11],gg[::11])]);want_a=np.array([m.gmul(int(g),int(a)^m.gmul(node,int(b)))^int(o) for g,a,b,o in zip(gg,aa,bb,A)]);want_b=np.array([m.gmul(int(g),int(b))^int(o) for g,b,o in zip(gg,bb,B)]);gout=np.array([m.gmul(m.ginv(int(g)),1 if (boundary and p==0) else m.gmul(int(v) or 1,int(v) or 1)) for g,v,p in zip(gg,want_b,pb)]);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(node,False,boundary);ci=build(node,True,boundary)
  for inv,c,src,dst in [(False,cf,[aa,bb,gg,pa,pb,ppa,ppb],[want_a,want_b,gout,pa,pb,ppa,ppb]),(True,ci,[want_a,want_b,gout,pa,pb,ppa,ppb],[aa,bb,gg,pa,pb,ppa,ppb])]:
   sim=km.Simulator(batch_size=shots,seed=142059+inv);sim.use_same_registers_as(c);names=['a','b','g','prev_a','prev_b','prior_a','prior_b']
   for name,v in zip(names,src):sim.write_across_shots(c.register_data[name],bits(v))
   sim.do(c);errors={}
   for name,v in zip(names,dst):
    bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(v)).any(0))
    if len(bad):errors[name]=bad[:8].tolist()
   public={q for w in c.register_data.values() for q in w};sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(boundary=boundary,inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
   sim.do(cf if inv else ci)
   for name,v in zip(names,src):assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),name
   assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',shots=shots*6,rows=rows,scope='arbitrary frame/pivot inputs, forced current/previous/prior zeros, three boundary charts, both directions/compositions, clean scratch and phase');Path('results/schur_gap_pivot_native.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':test()
