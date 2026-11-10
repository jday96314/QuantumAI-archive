"""Reversible completion-aware conversion from Schur heads to affine factors.

Retains the completed-prefix length, a boundary reciprocal, and an endpoint
ratio. These outputs are explicit; their width and cleanup are charged.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from base_inverse import inplace as invert12
from short_adder import add
from schur_control import caxpy


def build(N,inverse=False):
 m=setup()[0];m.inv_word=bind(m);n=N//2;b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(N)];degree=b.create_quantum_register(n.bit_length(),'degree');factor=b.create_quantum_register(12,'factor');seed=b.create_quantum_register(12,'seed');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 def count(back=False):
  for i in (reversed(range(n)) if back else range(n)):
   with m.eq_test(E,[q.id for q in S[2*i+1]],0) as z:
    E.x(z);add(b,km.array([km.q(z)]),degree,subtract=back);E.x(z)
 def selection():
  b.x(factor[0]);flag=b.alloc_qubits(1)[0]
  for i in range(n):
   b.flip_if_equal(degree,i+1,target=flag);caxpy(m,E,flag.id,[q.id for q in S[2*i+1]],[q.id for q in factor],1);b.cx(flag,factor[0]);b.flip_if_equal(degree,i+1,target=flag)
  b.free(km.array([flag]))
 def armor():
  flag=b.alloc_qubits(1)[0]
  for i in range(n):
   b.flip_if_less_than(degree,i+1,target=flag);b.cx(flag,S[2*i+1][0]);b.flip_if_less_than(degree,i+1,target=flag)
  b.free(km.array([flag]))
 def endpoint(back=False):
  with m.inv_word(E,[q.id for q in S[-1]]) as bi:
   rt=m.prod_rt([q.id for q in seed],[q.id for q in S[-2]],bi)
   if back:E.free(rt)
   else:E.mac(rt)
 def factors(back=False):
  for i in (range(n-1) if back else reversed(range(n-1))):
   olda,oldb=S[2*i:2*i+2];na,nb=S[2*i+2:2*i+4];B=[q.id for q in oldb]
   if not back:b.cx(olda,na)
   with m.inv_word(E,B) as bi,m.opened(E,'affine_parameter_product',12) as temp:
    for w in [na,nb]:m.mon_slot(E,[q.id for q in w],B if back else bi,bi if back else B,temp)
   if back:b.cx(olda,na)
 if inverse:factors(True);endpoint(True);armor();invert12(b,factor);selection();count(True)
 else:count();selection();invert12(b,factor);armor();endpoint();factors()
 ab.release_recycle('affine_parameter_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(N=6,shots=512):
 m=setup()[0];rng=np.random.default_rng(142065+N);n=N//2;source=[];want=[];dd=[];ff=[];rr=[]
 for case in range(shots):
  d=case%(n+1);v=rng.integers(0,4096,N).tolist();v[1::2]=rng.integers(1,4096,n).tolist();v[2*d:]=[0]*(N-2*d);src=v.copy();f=m.ginv(v[2*d-1]) if d else 1
  for i in range(d,n):v[2*i+1]=1
  r=m.gmul(v[-2],m.ginv(v[-1]))
  for i in reversed(range(n-1)):
   iv=m.ginv(v[2*i+1]);v[2*i+2]=m.gmul(v[2*i+2]^v[2*i],iv);v[2*i+3]=m.gmul(v[2*i+3],iv)
  source.append(src);want.append(v);dd.append(d);ff.append(f);rr.append(r)
 source=np.array(source);want=np.array(want);dd=np.array(dd);ff=np.array(ff);rr=np.array(rr);cf=build(N);ci=build(N,True);bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);rows=[]
 for inv,c,ss,tt in [(False,cf,source,want),(True,ci,want,source)]:
  sim=km.Simulator(batch_size=shots,seed=142065+inv);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(ss[:,i],12))
  if inv:
   for name,v in [('degree',dd),('factor',ff),('seed',rr)]:sim.write_across_shots(c.register_data[name],bits(v,len(c.register_data[name])))
  sim.do(c);errors={}
  for i in range(N):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(tt[:,i],12)).any(0))
   if len(bad):errors[f'S{i}']=bad[:8].tolist()
  for name,v in [('degree',dd),('factor',ff),('seed',rr)]:
   vv=np.zeros_like(v) if inv else v;bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(vv,len(c.register_data[name]))).any(0))
   if len(bad):errors[name]=bad[:8].tolist()
  public={q for w in c.register_data.values() for q in w};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(N=N,inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  sim.do(cf if inv else ci)
  for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(ss[:,i],12)),i
  assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
  c.to_file(f'candidates/schur_affine_parameters_n{N}_{"inverse" if inv else "forward"}.kmb',format='kmb')
 result=dict(status='PASS',N=N,shots=2*shots,rows=rows,scope='all prefix lengths, arbitrary nonzero active b and zero completed suffix, exact transformed fields/degree/boundary/seed, both directions/compositions, scratch and phase');Path(f'results/schur_affine_parameters_n{N}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.shots)
