"""Bijective two-word base Schur synthesis, including the completed zero pair."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from fast_dqi import FastBuilder
from schur_fold import axpy,xword


def build(kappa,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());a=b.create_quantum_register(12,'a');bv=b.create_quantum_register(12,'b');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));A=[q.id for q in a];B=[q.id for q in bv]
 def ratio(back=False):
  with m.eq_test(E,B,0) as z:
   E.cx(z,B[0])
   with m.inv_word(E,B) as bi,m.opened(E,'basefold_ratio_product',12) as temp:m.mon_slot(E,A,B if back else bi,bi if back else B,temp)
   E.cx(z,B[0])
 def normscale(back=False):
  with m.opened(E,'basefold_norm',12) as nn:
   nq=km.array([km.q(i) for i in nn]);cols=[m.gmul(1<<i,1<<i)^(1<<i) for i in range(12)]
   def linear():
    for i,q in enumerate(a):
     for j in range(12):
      if cols[i]>>j&1:b.cx(q,nq[j])
    xword(b,nq,kappa)
   linear()
   with m.inv_word(E,nn) as ni,m.opened(E,'basefold_norm_product',12) as temp:m.mon_slot(E,B,ni if back else nn,nn if back else ni,temp)
   linear()
 def shift():
  with m.eq_test(E,B,0) as z:b.x(a[0]);E.cx(z,A[0])
 if inverse:shift();normscale(True);ratio(True)
 else:ratio();normscale();shift()
 ab.release_recycle('basefold_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test():
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142047);shots=8192;aa=rng.integers(0,4096,shots);bb=rng.integers(0,4096,shots);bb[:4096]=0;aa[:4096]=np.arange(4096);r=np.array([m.gmul(int(a),m.ginv(int(b) or 1)) for a,b in zip(aa,bb)]);cc=np.array([m.gmul(int(b),m.gmul(int(t),int(t))^int(t)^k) for t,b in zip(r,bb)]);rr=r^(bb!=0);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(k);ci=build(k,True);rows=[]
 for rev,c,src,dst in [(False,cf,[aa,bb],[rr,cc]),(True,ci,[rr,cc],[aa,bb])]:
  sim=km.Simulator(batch_size=shots,seed=142047+rev);sim.use_same_registers_as(c)
  for name,v in zip(['a','b'],src):sim.write_across_shots(c.register_data[name],bits(v))
  sim.do(c)
  for name,v in zip(['a','b'],dst):assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),name
  public={q for word in c.register_data.values() for q in word};sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any();sim.do(cf if rev else ci)
  for name,v in zip(['a','b'],src):assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),name
  assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any();rows.append(dict(inverse=rev,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
 out=dict(status='PASS',shots=2*shots,rows=rows,scope='all 4096 zero-divisor inputs, random nonzero pivots, both directions/compositions, clean scratch and phase');Path('results/schur_basefold_native.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':test()
