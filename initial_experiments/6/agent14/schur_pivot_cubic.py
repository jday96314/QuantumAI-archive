"""Schur pivot conversion with paired scale and reciprocal-frame erasures."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from cubic_low_workspace_inverse import bind
from tower_inverse import apply_map,map_ops
from schur_affine_outer_seedless import ratio_phase
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create
from spooky_extension_norm_inverse import init_delayed
from cubic_mask_phase import native_phase_inverse as inverse_mask_native
from recursive_forney import inverse_dot_data


def transpose_native(m,b,word,mask,target):
 flag=b.alloc_bits(1)[0];b.bit_store(flag,False)
 for ll,rr,oo in m.KT:
  lead=km.q(word[ll[0]]);other=km.array([km.q(word[i]) for i in ll[1:]]);right=km.array([km.q(target[i]) for i in rr]);mm=km.array([mask[i] for i in oo])
  b.cx(other,lead);b.cx(mm,flag);b.push_condition(flag);b.cx(lead,right);b.pop_condition();b.cx(mm,flag);b.cx(other,lead)
 b.free(km.array([flag]))


def dot_inverse_native(m,b,source,u,fw,bk):
 source=km.array([km.q(q) for q in source]);u=km.array([km.q(q) for q in u]);sops=map_ops(bk);trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));uops=map_ops(trans)
 apply_map(b,source,sops);apply_map(b,u,uops);leaf,iv,cleanup=inverse_dot_data(b,source,u);b.cz(leaf,iv);cleanup();apply_map(b,u,uops,True);apply_map(b,source,sops,True)


def build(node,inverse=False,first=False):
 m,k,fw,bk,tr=setup();m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());lo=b.create_quantum_register(12,'a');hi=b.create_quantum_register(12,'b');g=b.create_quantum_register(12,'g');pa=b.create_quantum_register(12,'prev_a');pb=b.create_quantum_register(12,'prev_b');ghost=b.create_classical_register(36,'ghost');ab=Adapter(b);sc,E=lean_create(m,ab);A=[q.id for q in lo];B=[q.id for q in hi];G=[q.id for q in g]
 def node_axpy():
  for i,q in enumerate(hi):
   c=m.gmul(node,1<<i)
   for j in range(12):
    if c>>j&1:b.cx(q,lo[j])
 def offset():
  if not first:b.cx(pa,lo);b.cx(pb,lo);b.cx(pb,hi)
 def scale(back=False):
  if back:
   mask=b.alloc_bits(24)
   with m.opened(E,'pivot_scale_values',12) as tmp:
    for p,word in enumerate([A,B]):
     E.mac(m.prod_rt(tmp,word,G));b.hmr(km.array([km.q(q) for q in word]),mask[12*p:12*p+12])
     for q,r in zip(word,tmp):E.b.swap(q,r)
     for q,bit in zip(word,ghost[12*p:12*p+12]):ab.z(q,bit)
   with m.opened(E,'pivot_scale_transpose',12) as u:
    for p,word in enumerate([A,B]):transpose_native(m,b,word,mask[12*p:12*p+12],u)
    dot_inverse_native(m,b,G,u,fw,bk)
    for p,word in reversed(list(enumerate([A,B]))):transpose_native(m,b,word,mask[12*p:12*p+12],u)
   b.bit_store(mask,False);b.free(mask)
  else:
   with m.opened(E,'pivot_scale_values',12) as out:
    for p,word in enumerate([A,B]):
     E.mac(m.prod_rt(out,G,word));b.hmr(km.array([km.q(q) for q in word]),ghost[12*p:12*p+12])
     for q,r in zip(word,out):E.b.swap(q,r)
 def gupdate():
  z=b.alloc_qubits(1)[0].id
  with m.eq_test(E,B,0) as zz:E.cx(zz,z)
  if True:
   E.cx(z,B[0]);ops=map_ops([m.gmul(1<<i,1<<i) for i in range(12)]);apply_map(b,hi,ops)
   if inverse:
    with m.opened(E,'pivot_old_regular_frame',12) as prod:
     rt=m.prod_rt(prod,B,G);E.mac(rt);out=init_delayed(m,E,prod)
     for q,bit in zip(prod,ghost[24:]):ab.z(q,bit)
     E.free(rt)
    mask=b.alloc_bits(12);b.hmr(g,mask)
    for q,r in zip(G,out):E.b.swap(q,r)
    E.A.close('ext_delayed_inverse',out)
    with m.opened(E,'pivot_reciprocal_repay',12) as prod:
     rt=m.prod_rt(prod,B,G);E.mac(rt);inverse_mask_native(b,prod,mask,fw,bk);E.free(rt)
    b.bit_store(mask,False);b.free(mask)
   else:
    with m.opened(E,'pivot_new_frame',12) as out:
     with m.inv_word(E,G) as gi:E.mac(m.prod_rt(out,B,gi))
     b.hmr(g,ghost[24:])
     for q,r in zip(G,out):E.b.swap(q,r)
   apply_map(b,hi,ops,True);E.cx(z,B[0])
  with m.eq_test(E,B,0) as zz:E.cx(zz,z)
  b.free(km.array([km.q(z)]))
 if inverse:gupdate();offset();scale(True);node_axpy();b.bit_store(ghost,False)
 else:b.bit_store(ghost,False);node_axpy();scale();offset();gupdate()
 ab.release_recycle('pivot_spooky_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(36,c.num_bits)))


def test():
 m,k,fw,bk,tr=setup();node=tr(26,fw);rng=np.random.default_rng(142097);shots=4096;aa,bb,pa,pb=rng.integers(0,4096,(4,shots));gg=rng.integers(1,4096,shots);bb[::7]=np.array([m.gmul(int(p),m.ginv(int(g))) for p,g in zip(pb[::7],gg[::7])]);rawa=np.array([int(a)^m.gmul(node,int(v)) for a,v in zip(aa,bb)]);a=np.array([m.gmul(int(g),int(v))^int(p)^int(q) for g,v,p,q in zip(gg,rawa,pa,pb)]);v=np.array([m.gmul(int(g),int(x))^int(p) for g,x,p in zip(gg,bb,pb)]);gout=np.array([m.gmul(m.gmul(int(x) or 1,int(x) or 1),m.ginv(int(g))) for g,x in zip(gg,v)]);bits=lambda x:((x[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(node);ci=build(node,True);rows=[]
 h=np.array([m.ginv(int(x)) for x in gg]);hout=np.array([m.ginv(int(x)) for x in gout])
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142097+inv);sim.use_same_registers_as(c);source=[a,v,hout,pa,pb] if inv else [aa,bb,gg,pa,pb];want=[aa,bb,h,pa,pb] if inv else [a,v,gout,pa,pb];names=['a','b','g','prev_a','prev_b']
  for name,x in zip(names,source):sim.write_across_shots(c.register_data[name],bits(x))
  sim.do(c);errors=[name for name,x in zip(names,want) if not np.array_equal(sim.read_across_shots(c.register_data[name]),bits(x))];mask=sim.read_across_shots(c.register_data['ghost']);phase=np.zeros(shots,dtype=bool) if inv else ((mask&np.concatenate([bits(rawa),bits(bb),bits(gg)])).sum(0)%2).astype(bool);public={q for name,reg in c.register_data.items() if name!='ghost' for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(scratch).any()),phase=bool((sim.read_phase_flipped_across_shots()!=phase).any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
  if not inv:
   sim.write_across_shots(c.register_data['g'],bits(hout))
   sim.do(ci)
   sim.write_across_shots(c.register_data['g'],bits(gg))
   for name,x in zip(names,source):assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(x))
   assert not sim.read_across_shots(scratch).any();assert not sim.read_across_shots(c.register_data['ghost']).any();assert not sim.read_phase_flipped_across_shots().any()
 Path('results/schur_pivot_cubic_native.json').write_text(json.dumps(dict(status='PASS',shots=shots*2,rows=rows,scope='arbitrary nonzero frame, zero-armed current b; reverse frame reciprocal; exact held phases and fresh inverse. Pair host changes frame basis between stages, full synthesis bridge still required.'),indent=2)+'\n')
if __name__=='__main__':test()
