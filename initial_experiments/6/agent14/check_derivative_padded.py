"""Native consistent singular solves on every physical padding length."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from derivative_inverse_padded import build
from extension_inverse import setup
from quadratic_schur_model import F


def test(n,shots=512):
 m,k,fw,bk,tr=setup();rng=np.random.default_rng(142098+n);cases=[]
 for row in range(shots):
  d=row%(n+1);roots=rng.choice(np.arange(1,4096),d,replace=False).tolist();B=[1]
  for a in roots:B=F.product(B,[a,1])
  A=rng.integers(0,4096,d).tolist();c=n-d;cases.append(dict(A=[0]*c+A,B=[0]*c+B[:-1],roots=roots,values=[F.mul(F.evaluate(A,a),F.inv(F.evaluate(F.derivative(B),a))) for a in roots],c=c))
 bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool);cf=build(n);ci=build(n,True);sim=km.Simulator(batch_size=shots,seed=142098);sim.use_same_registers_as(cf)
 for i in range(n):
  for key in ['A','B']:sim.write_across_shots(cf.register_data[f'{key}{i}'],bits([tr(v[key][i],fw) for v in cases]))
 padding=bits([v['c']+1 for v in cases],len(cf.register_data['padding']));sim.write_across_shots(cf.register_data['padding'],padding);sim.do(cf);out=np.array([sim.read_across_shots(cf.register_data[f'A{i}']).astype(np.int64).T@(1<<np.arange(12)) for i in range(n)]).T;errors=[]
 for row,(v,w) in enumerate(zip(cases,out)):
  W=[tr(int(x),bk) for x in w]
  if any(F.evaluate(W,a)!=want for a,want in zip(v['roots'],v['values'])):errors.append(row)
 public={q for reg in cf.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(max(cf.num_qubits,ci.num_qubits)) if km.q(i) not in public]);r=dict(n=n,Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like,errors=errors[:12],scratch=bool(sim.read_across_shots(scratch).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));print(json.dumps(r),flush=True)
 assert not errors and not r['scratch'] and not r['phase'];assert np.array_equal(sim.read_across_shots(cf.register_data['padding']),padding)
 for i in range(n):assert np.array_equal(sim.read_across_shots(cf.register_data[f'B{i}']),bits([tr(v['B'][i],fw) for v in cases]))
 sim.do(ci)
 for i in range(n):
  for key in ['A','B']:assert np.array_equal(sim.read_across_shots(ci.register_data[f'{key}{i}']),bits([tr(v[key][i],fw) for v in cases]))
 assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any();assert np.array_equal(sim.read_across_shots(ci.register_data['padding']),padding);r.update(status='PASS',shots=shots,scope='all common monomial padding lengths, arbitrary consistent numerator, physical residue values, restored B/count/scratch/phase and exact inverse composition; rare leading-gap and evaluation-at-one poles charged separately');Path(f'results/derivative_inverse_padded_n{n}_native.json').write_text(json.dumps(r,indent=2)+'\n')
if __name__=='__main__':
 for n in [3,7,35]:test(n)
