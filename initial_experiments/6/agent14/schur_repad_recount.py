"""Reconstruct common-padding multiplicity from known irreducible factors.

This applies when the true denominator splits over the nonzero base field,
so none of the fixed quadratic padding factors occurs in the true part.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from schur_dynamic_fold import powers
from schur_fold import axpy,xword
from quadratic_schur_model import F,node_maps
from fast_dqi import FastBuilder


def build(n,kind='recount'):
 m,k,fw,bk,tr=setup();b=FastBuilder(km.CircuitBuilder());B=[b.create_quantum_register(12,f'B{i}') for i in range(n)];index=b.create_quantum_register((n+1).bit_length(),'index');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));ix=[q.id for q in index]
 if kind!='recount':
  if kind=='decrement':b.x(index)
  with m.opened(E,'recount_increment',len(ix)-2) as tree:E.inc1(ix,tree)
  if kind=='decrement':b.x(index)
 else:
  _,kap,_,_=node_maps(n)
  with m.eq_test(E,[q.id for q in B[0]],0) as z:E.cx(z,ix[0])
  for threshold in range(2,n,2):
   c=tr(kap[(n-threshold-1)//2],fw);pp=powers(n,c,m)
   def reduce():
    for j in range(2,n):axpy(b,B[j],B[0],pp[j][0],m);axpy(b,B[j],B[1],pp[j][1],m)
    xword(b,B[0],pp[n][0]);xword(b,B[1],pp[n][1])
   reduce()
   with m.eq_test(E,[q.id for word in B[:2] for q in word],0) as z:
    for bit,q in enumerate(ix):
     if (threshold^(threshold-2))>>bit&1:E.cx(z,q)
   reduce()
 ab.release_recycle('recount_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(n,shots=512):
 from schur_repad import build as padding
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142096+n);aa=[];bb=[];counts=[];_,kap,_,_=node_maps(n)
 for case in range(shots):
  d=case%(n+1);roots=rng.choice(np.arange(1,4096),d,replace=False);B=[1]
  for root in roots:B=F.product(B,[int(root),1])
  c=n-d;P=[0,1] if c%2 else [1]
  for kk in range(1,n,2):
   if kk>=d:P=F.product(P,[kap[(kk-1)//2],1,1])
  bb.append(F.product(P,B)[:-1]);counts.append(c)
 bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool);c=build(n);sim=km.Simulator(batch_size=shots,seed=142096);sim.use_same_registers_as(c)
 for i in range(n):sim.write_across_shots(c.register_data[f'B{i}'],bits([tr(v[i],fw) for v in bb]))
 sim.do(c);assert np.array_equal(sim.read_across_shots(c.register_data['index']),bits(counts,len(c.register_data['index'])))
 for i in range(n):assert np.array_equal(sim.read_across_shots(c.register_data[f'B{i}']),bits([tr(v[i],fw) for v in bb]))
 public={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any();sim.do(c);assert not sim.read_across_shots(c.register_data['index']).any();assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any();r=dict(status='PASS',n=n,shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,scope='all padding lengths, random distinct nonzero physical roots, exact reconstructed counter, restored B, double composition and clean scratch/phase');Path(f'results/schur_repad_recount_n{n}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 for n in [3,7,35]:test(n)
