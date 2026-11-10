"""Recreate the Schur affine frame from the original field pivots."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup,Adapter
from low_workspace_inverse import bind
from tower_inverse import map_ops,apply_map
from fast_dqi import FastBuilder


def build(N):
 m=setup()[0];m.inv_word=bind(m);n=N//2;b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(N)];g=b.create_quantum_register(12,'g');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));G=[q.id for q in g];ops=map_ops([m.gmul(1<<i,1<<i) for i in range(12)]);b.x(g[0])
 for i in range(n):
  word=S[2*i+1];B=[q.id for q in word]
  with m.eq_test(E,B,0) as z:
   E.cx(z,B[0])
   with m.inv_word(E,B) as bi,m.opened(E,'frame_product',12) as temp:
    ib=km.array([km.q(q) for q in bi]);apply_map(b,word,ops);apply_map(b,ib,ops)
    back=(n-1-i)%2!=0;m.mon_slot(E,G,bi if back else B,B if back else bi,temp)
    apply_map(b,ib,ops,True);apply_map(b,word,ops,True)
   E.cx(z,B[0])
 ab.release_recycle('frame_product_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(N):
 m=setup()[0];shots=512;n=N//2;rng=np.random.default_rng(142066+N);ss=rng.integers(0,4096,(shots,N));ss[:,1::2]=rng.integers(1,4096,(shots,n));want=[]
 for j in range(shots):
  d=j%(n+1);ss[j,2*d:]=0;g=1
  for i in range(n):g=m.gmul(m.gmul(int(ss[j,2*i+1]) or 1,int(ss[j,2*i+1]) or 1),m.ginv(g))
  want.append(g)
 c=build(N);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142066);sim.use_same_registers_as(c)
 for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(ss[:,i]))
 sim.do(c);assert np.array_equal(sim.read_across_shots(c.register_data['g']),bits(want))
 for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(ss[:,i]))
 public={q for w in c.register_data.values() for q in w};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any();c.to_file(f'candidates/schur_frame_product_n{N}.kmb',format='kmb');r=dict(status='PASS',N=N,shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,scope='all completed-prefix lengths and arbitrary active fields, exact frame recurrence, restored input and clean scratch/phase');Path(f'results/schur_frame_product_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
