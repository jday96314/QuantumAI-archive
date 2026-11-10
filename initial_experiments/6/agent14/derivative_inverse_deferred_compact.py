"""Defer companion scalar measurements until the matching inverse transform."""
import argparse,json
from pathlib import Path
import kickmix as km
from extension_inverse import Adapter,setup
from cubic_low_workspace_inverse_compact import bind
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create
from inverse_ratio_mask_phase import apply as ratio_phase
import derivative_inverse_raw as raw
import bezout_payload_monic as monic
import bezout_payload_batched as batched


def build(n,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());aa=[b.create_quantum_register(12,f'A{i}') for i in range(n)];bb=[b.create_quantum_register(12,f'B{i}') for i in range(n)];ghost=b.create_classical_register(12*(n-1),'companion_ghost');ab=Adapter(b);sc,E=lean_create(m,ab)
 def companion(m,E,f,payloads,flag,inverse=False,one=False):
  size=len(payloads[0]);offset=0 if one else 24*(n//2-size+1);slots=[ghost[offset+12*j:offset+12*(j+1)] for j in range(2)]
  def rotate(back=False):
   for y in payloads:
    for i in (range(size-1) if back else reversed(range(size-1))):
     for q,r in zip(y[i],y[i+1]):b.cswap(km.q(flag),km.q(q),km.q(r))
  if not inverse:rotate()
  source=f[0] if inverse else ([] if one else f[-1])
  from contextlib import nullcontext
  with (m.inv_word(E,source) if source else nullcontext(None)) as si:
   for y,held in zip(payloads,slots):
    masked=E.A.open('deferred_companion_mask',12)
    if True:
     E.mac([('P',(flag,),0,(q,),0,(r,)) for q,r in zip(y[0],masked)])
     if si is not None:
      with m.opened(E,'deferred_companion_scale',12) as tmp:m.mon_slot(E,masked,si,source,tmp)
     if inverse:
      for q,bit in zip(masked,held):ab.z(q,bit)
      b.bit_store(held,False)
     rt=[];m.bil_terms(rt,m.KT,m.spec_plain(masked),m.spec_plain(f[0]) if one else [a^c for a,c in zip(m.spec_plain(f[0]),m.spec_plain(f[-1]))],y[0]);E.mac(rt)
     if one:
      for q,r in zip(masked,y[0]):E.cx(q,r)
     for i in range(1,size):E.mac(m.prod_rt(y[i],masked,f[i]))
     if inverse:
      fresh=b.alloc_bits(12);b.hmr(km.array([km.q(q) for q in masked]),fresh);E.A.close('deferred_companion_mask',masked)
      if one:
       for q,bit in zip(y[0],fresh):b.push_condition(bit);b.cz(km.q(flag),km.q(q));b.pop_condition()
      else:ratio_phase(b,f[-1],y[0],flag,fresh)
      b.bit_store(fresh,False);b.free(fresh)
     else:
      b.hmr(km.array([km.q(q) for q in masked]),held);E.A.close('deferred_companion_mask',masked)
  if inverse:rotate(True)
 oldm,oldb=monic.companion,batched.companion
 try:
  monic.companion=lambda m,E,f,payloads,flag,inverse=False:companion(m,E,f,payloads,flag,inverse,True)
  batched.companion=companion
  raw.apply(m,E,[[q.id for q in w] for w in aa],[[q.id for q in w] for w in bb],[],inverse)
 finally:monic.companion,batched.companion=oldm,oldb
 ab.release_recycle('deferred_emitter',sc);c=b.finish_circuit();held={bit.id for bit in ghost};return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits) if i not in held))


def test(n,shots=256):
 # The established independent polynomial model remains unchanged.  Its
 # forward/pair test below additionally accounts for the held classical masks.
 import numpy as np
 from derivative_inverse_model import invmod,divide
 from quadratic_schur_model import F
 m,_,fw,bk,tr=setup();rng=np.random.default_rng(142122+n);cases=[]
 while len(cases)<shots:
  roots=rng.choice(np.arange(1,4096),n,replace=False).tolist();B=[1]
  for a in roots:B=F.product(B,[int(a),1])
  f=[F.pow(x,2048) for x in B[::2]];g=[F.pow(x,2048) for x in B[1::2]]
  if len(cases)%2 and n>=7:
   f=rng.integers(0,4096,n//2+1).tolist();g=rng.integers(0,4096,n//2).tolist()+[1];f[-2]=F.mul(f[-1],g[-2]);B=[0]*(n+1)
   for i,v in enumerate(f):B[2*i]=F.mul(v,v)
   for i,v in enumerate(g):B[2*i+1]=F.mul(v,v)
  _,r=divide(f,g);f,g=g,r;good=True;used=False;expected=n//2
  while expected:
   if len(g)!=expected:
    if not used and len(g)==expected-1 and f[0]:g=[0]+g;used=True
    else:good=False;break
   if expected==1:break
   _,r=divide(f,g);f,g=g,r;expected-=1
  if not good or not g or not g[0]:continue
  A=rng.integers(0,4096,n).tolist();want=divide(F.product(A,invmod(F.derivative(B),B)),B)[1];want+=[0]*(n-len(want));cases.append((A,B[:-1],want))
 bits=lambda v:((np.array(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(n);ci=build(n,True);sim=km.Simulator(batch_size=shots,seed=142122);sim.use_same_registers_as(cf)
 for k,bank in [('A',0),('B',1)]:
  for i in range(n):sim.write_across_shots(cf.register_data[f'{k}{i}'],bits([tr(v[bank][i],fw) for v in cases]))
 sim.do(cf)
 for k,bank in [('A',2),('B',1)]:
  for i in range(n):assert np.array_equal(sim.read_across_shots(cf.register_data[f'{k}{i}']),bits([tr(v[bank][i],fw) for v in cases])),(k,i)
 sim.do(ci)
 for k,bank in [('A',0),('B',1)]:
  for i in range(n):assert np.array_equal(sim.read_across_shots(cf.register_data[f'{k}{i}']),bits([tr(v[bank][i],fw) for v in cases])),(k,i)
 public={q for name,reg in cf.register_data.items() if name!='companion_ghost' for q in reg};junk=km.array([km.q(q) for q in range(max(cf.num_qubits,ci.num_qubits)) if km.q(q) not in public]);assert not sim.read_across_shots(junk).any();assert not sim.read_across_shots(cf.register_data['companion_ghost']).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS',n=n,shots=shots,forward=dict(Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like),inverse=dict(Q=ci.num_qubits,T=ci.operation_statistics().num_toffoli_like));Path(f'results/derivative_inverse_deferred_compact_n{n}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=7);p.add_argument('--shots',type=int,default=256);a=p.parse_args();test(a.n,a.shots)
