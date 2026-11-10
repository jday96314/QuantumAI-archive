"""Native comparison of deferred head erasure against the established step."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from schur_affine_dynamic_step_base import build
from schur_affine_dynamic_step_canonical import build as old_build
from extension_inverse import setup
from quadratic_schur_model import F


def test(d,shots=512):
 m,k,fw,bk,tr=setup();deltas=[tr(6^F.mul(6,6),fw),tr(14^F.mul(14,14),fw)];old=old_build(d,k,deltas);cf=build(d,k,deltas);ci=build(d,k,deltas,True);rng=np.random.default_rng(142090+d);bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool)
 src={name:rng.integers(0,4096,shots) for name in old.register_data if name!='flags'};src['t']=rng.integers(1,4096,shots);src['Nlead']=np.zeros(shots,dtype=int);src['Btop']=np.zeros(shots,dtype=int);comp=np.arange(shots)%3==0
 for name in src:
  if name.startswith('B') or (name.startswith('N') and name!=f'N{d}'):src[name][comp]=0
 src['u'][comp]=0;src['t'][comp]=1;choice=rng.integers(0,3,shots);flags=np.array([choice==1,choice==2]);flags[:,comp]=False
 baseline=km.Simulator(batch_size=shots,seed=142090);baseline.use_same_registers_as(old)
 for name,v in src.items():baseline.write_across_shots(old.register_data[name],bits(v))
 baseline.write_across_shots(old.register_data['flags'],flags);baseline.do(old);assert not baseline.read_phase_flipped_across_shots().any();expected={name:baseline.read_across_shots(reg) for name,reg in old.register_data.items()}
 rows=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142091+inv);sim.use_same_registers_as(c)
  for name,v in src.items():sim.write_across_shots(c.register_data[name],expected[name] if inv else bits(v))
  sim.write_across_shots(c.register_data['flags'],flags);sim.do(c);errors=[]
  for name in src:
   if not np.array_equal(sim.read_across_shots(c.register_data[name]),bits(src[name]) if inv else expected[name]):errors.append(name)
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),flags)
  public={q for name,reg in c.register_data.items() if name!='head_ghost' for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);ghost=sim.read_across_shots(c.register_data['head_ghost']);wantphase=np.zeros(shots,dtype=bool) if inv else ((ghost[:24]&np.concatenate([bits(src['u']),bits(src['t'])])).sum(0)%2!=0)^(ghost[24]&comp)
  r=dict(d=d,inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(scratch).any()),phase=bool((sim.read_phase_flipped_across_shots()!=wantphase).any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
  if not inv:
   sim.do(ci)
   for name in src:assert np.array_equal(sim.read_across_shots(ci.register_data[name]),bits(src[name]))
   assert not sim.read_across_shots(scratch).any();assert not sim.read_across_shots(ci.register_data['head_ghost']).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',shots=shots*2,rows=rows,scope='arbitrary generic and completed tails, mutually exclusive node flags, exact old-step output comparison, forward held head/completion phase, fresh inverse and paired cleanup');Path(f'results/schur_affine_dynamic_step_base_d{d}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 for d in [0,3,33]:test(d)
