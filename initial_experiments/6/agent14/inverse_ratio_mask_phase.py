"""Mask-selected parity of flag * numerator / divisor, without a quotient word.

A runtime classical mask scales the divisor through reversible row operations.
The tower inverse-dot circuit then gates its final three-bit pairing by flag.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from cubic_mask_phase_joint import classical_scale_inverse
from cubic_trace_field import data
from extension_inverse import setup
from recursive_forney import inverse_dot_data
from tower_inverse import map_ops,apply_map
from fast_dqi import FastBuilder


def apply(b,divisor,numerator,flag,mask):
 m,k,fw,bk,tr=setup();field=km.GF2Field(12)
 source=km.array([km.q(q) for q in divisor]);num=km.array([km.q(q) for q in numerator]);bops=map_ops(bk)
 trans=tuple(sum(((fw[i]>>j)&1)<<i for i in range(12)) for j in range(12));mops=map_ops(trans)
 def trace(v):
  ans=0
  for _ in range(12):ans^=v;v=field.mul(v,v)
  assert ans in (0,1);return ans
 gram=map_ops([sum(trace(field.mul(1<<i,1<<j))<<i for i in range(12)) for j in range(12)])
 active=b.alloc_bits(1)[0];b.bit_store(active,False);b.bit_store(active,True,control=mask);b.push_condition(active)
 apply_map(b,source,bops);apply_map(b,num,bops);apply_map(b,mask,mops)
 dual=map_ops(data()['dual_columns']);apply_map(b,mask,dual);record=classical_scale_inverse(b,source,mask);apply_map(b,num,gram)
 leaf,iv,cleanup=inverse_dot_data(b,source,num)
 for q,r in zip(leaf,iv):b.ccz(km.q(flag),q,r)
 cleanup();apply_map(b,num,gram,True)
 for kind,j,i,bit in reversed(record):
  b.push_condition(bit);getattr(b,kind)(source[j],source[i]);b.pop_condition();b.bit_store(bit,False);b.free(km.array([bit]))
 apply_map(b,mask,dual,True);apply_map(b,mask,mops,True);apply_map(b,num,bops,True);apply_map(b,source,bops,True)
 b.pop_condition();b.bit_store(active,False);b.free(km.array([active]))


def test():
 b=FastBuilder(km.CircuitBuilder());d=b.create_quantum_register(12,'divisor');n=b.create_quantum_register(12,'numerator');flag=b.create_quantum_register(1,'flag');mask=b.create_classical_register(12,'mask');apply(b,[q.id for q in d],[q.id for q in n],flag[0].id,mask);c=b.finish_circuit();m=setup()[0];rng=np.random.default_rng(142121);shots=16384;dd=np.tile(np.arange(4096),4);nn=rng.integers(0,4096,shots);mm=rng.integers(0,4096,shots);ff=np.repeat([0,1],8192);mm[::11]=0;nn[::13]=0;want=np.array([bool(f and ((m.gmul(int(y),m.ginv(int(x)))&int(z)).bit_count()%2)) for x,y,z,f in zip(dd,nn,mm,ff)]);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142121);sim.use_same_registers_as(c)
 for reg,v in [(d,dd),(n,nn),(mask,mm),(flag,ff)]:sim.write_across_shots(reg,bits(v,len(reg)))
 sim.do(c)
 for reg,v in [(d,dd),(n,nn),(mask,mm),(flag,ff)]:assert np.array_equal(sim.read_across_shots(reg),bits(v,len(reg)))
 public=set(d)|set(n)|set(flag);junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);bad=np.flatnonzero(sim.read_phase_flipped_across_shots()!=want);r=dict(status='PASS' if not len(bad) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=shots,phase_errors=bad[:10].tolist(),scratch=bool(sim.read_across_shots(junk).any()));print(json.dumps(r),flush=True);Path('results/inverse_ratio_mask_phase_native.json').write_text(json.dumps(r,indent=2)+'\n');assert r['status']=='PASS' and not r['scratch']
if __name__=='__main__':test()
