"""Prepare a projective affine terminal tail from the infinity jet."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from lean_emitter import create as lean_create


def build(N,inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());heads=[b.create_quantum_register(12,f'S{i}') for i in range(N)];ghost=b.create_classical_register(49,'endpoint_ghost');ab=Adapter(b);sc,E=lean_create(m,ab);jet=[b.alloc_qubits(12) for _ in range(4)];p0,p1,q0,q1=[[q.id for q in w] for w in jet];flag=b.alloc_qubits(1)[0];n=N//2
 def step(i,back=False):
  aa=[q.id for q in heads[2*i]];bb=[q.id for q in heads[2*i+1]]
  def update():E.mac(m.prod_rt(p0,bb,q0));E.mac(m.prod_rt(p1,bb,q1));E.mac(m.prod_rt(p1,aa,q0))
  def swap():
   for a,c in [(p0,q0),(p1,q1)]:
    for x,y in zip(a,c):ab.swap(x,y)
  def lin():
   for x,y in zip(p0,p1):ab.cx(x,y)
  if back:lin();swap();update()
  else:update();swap();lin()
 def run(back=False):
  if not back:ab.x(q0[0])
  for i in (reversed(range(n)) if back else range(n)):step(i,back)
  if back:ab.x(q0[0])
 def predicate():
  with m.eq_test(E,q0,0) as z:ab.cx(z,flag.id)
 if not inverse:run();predicate()
 else:
  # Match the forward predicate's released scratch order before assigning
  # the public tail ports, so both directions use identical physical ports.
  with m.opened(E,'endpoint_port_reservation',11):pass
 tail=[b.alloc_qubits(12) for key in range(3)];headers=[]
 for i,(key,word) in enumerate(zip(['tailN0','tailN1','tailB0'],tail)):
  headers.extend(f'APPEND_TO_REGISTER q{q.id} r{N+1+i}' for q in word);headers.append(f'REGISTER r{N+1+i} \"{key}\"')
 t0,t1,d0=[[q.id for q in w] for w in tail];rt0=m.prod_rt(t0,p1,q0)+m.prod_rt(t0,p0,q1)+[('L',(flag.id,),0,(t0[0],))];rt1=m.prod_rt(t1,p0,q0);rtd=[];m.lin_rows_terms(rtd,[sum(((m.gmul(1<<i,1<<i)>>j)&1)<<i for i in range(12)) for j in range(12)],q0,d0)
 if inverse:
  run()
  for w,gg in zip(jet,[ghost[12*i:12*i+12] for i in range(4)]):
   for q,bit in zip(w,gg):ab.z(q.id,bit)
  E.free(rt1);E.free(rtd);b.free(tail[1]);b.free(tail[2]);predicate();ab.z(flag.id,ghost[48]);E.free(rt0);b.free(tail[0]);predicate();b.free(km.array([flag]));run(True)
  for w in jet:b.free(w)
  b.bit_store(ghost,False)
 else:
  E.mac(rt0);E.mac(rt1);E.mac(rtd);b.hmr(flag,ghost[48]);b.free(km.array([flag]));b.hmr(km.array([q for w in jet for q in w]),ghost[:48])
  for w in jet:b.free(w)
 ab.release_recycle('raw_endpoint_emitter',sc);c=b.finish_circuit();return km.Circuit('\n'.join(headers)+'\n'+str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(49,c.num_bits)))


def test(N=6,shots=512):
 from raw_cf_model import extract
 from quadratic_schur_model import F
 m,_,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142115+N);src=rng.integers(0,4096,(shots,N));src[0]=0
 for case in range(0,shots,2):
  d=(case//2)%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();hh,_=extract(F.syndrome(pos,val,N));src[case]=[tr(x,fw) for h in hh for x in h]
 jets=[];tails=[]
 for row in src:
  p0=p1=q1=0;q0=1
  for a,h in row.reshape(n,2):p0,p1,q0,q1=q0,q0^q1,p0^m.gmul(int(h),q0),p1^m.gmul(int(h),q1)^m.gmul(int(a),q0)
  jets.append([p0,p1,q0,q1]);tails.append([m.gmul(p1,q0)^m.gmul(p0,q1)^int(not q0),m.gmul(p0,q0),m.gmul(q0,q0)])
 jets=np.array(jets);tails=np.array(tails);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(N);ci=build(N,True);rows=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142115+inv);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  if inv:
   for i,key in enumerate(['tailN0','tailN1','tailB0']):sim.write_across_shots(c.register_data[key],bits(tails[:,i]))
  sim.do(c);errors=[]
  for i in range(N):
   if not np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i])):errors.append(f'S{i}')
  if not inv:
   for i,key in enumerate(['tailN0','tailN1','tailB0']):
    if not np.array_equal(sim.read_across_shots(c.register_data[key]),bits(tails[:,i])):errors.append(key)
  used={q for name,reg in c.register_data.items() if name!='endpoint_ghost' and (not inv or not name.startswith('tail')) for q in reg};scratch=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in used]);mask=sim.read_across_shots(c.register_data['endpoint_ghost']);expected=np.zeros(shots,dtype=bool) if inv else ((mask[:48]&np.concatenate([bits(jets[:,i]) for i in range(4)])).sum(0)%2).astype(bool)^(mask[48]&(jets[:,2]==0));r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(scratch).any()),phase=bool((sim.read_phase_flipped_across_shots()!=expected).any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
  if not inv:
   sim.do(ci)
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i]))
   sc=km.array([km.q(q) for q in range(max(cf.num_qubits,ci.num_qubits)) if q>=12*N]);assert not sim.read_across_shots(sc).any();assert not sim.read_across_shots(c.register_data['endpoint_ghost']).any();assert not sim.read_phase_flipped_across_shots().any()
 Path(f'results/raw_cf_endpoint_n{N}_native.json').write_text(json.dumps(dict(status='PASS',shots=shots,rows=rows),indent=2)+'\n')
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
