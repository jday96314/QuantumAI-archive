"""Reversible three-field infinity jet, without affine reciprocal steps."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from lean_emitter import create
from tower_inverse import map_ops,apply_map


def build(inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());regs={name:b.create_quantum_register(12,name) for name in ['a','b','p','q','d']};ab=Adapter(b);sc,E=create(m,ab);w={k:[q.id for q in v] for k,v in regs.items()};sq=map_ops([m.gmul(1<<i,1<<i) for i in range(12)])
 def swap():b.swap(regs['p'],regs['q'])
 def top():E.mac(m.prod_rt(w['p'],w['b'],w['q']))
 def deriv():
  E.mac(m.prod_rt(w['d'],w['q'],w['p']));apply_map(b,regs['q'],sq);E.mac(m.prod_rt(w['d'],w['a'],w['q']));apply_map(b,regs['q'],sq,True)
 if inverse:swap();deriv();top()
 else:top();deriv();swap()
 ab.release_recycle('three_jet_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(shots=8192):
 m=setup()[0];rng=np.random.default_rng(142137);aa,bb,pp,qq,dd=rng.integers(0,4096,(5,shots));pp[:4096]=np.arange(4096);qq[::7]=0;wantp=qq;wantq=np.array([int(p)^m.gmul(int(b),int(q)) for p,b,q in zip(pp,bb,qq)]);wantd=np.array([int(d)^m.gmul(int(p),int(q))^m.gmul(int(a)^int(b),m.gmul(int(q),int(q))) for a,b,p,q,d in zip(aa,bb,pp,qq,dd)]);src=dict(a=aa,b=bb,p=pp,q=qq,d=dd);dst=dict(a=aa,b=bb,p=wantp,q=wantq,d=wantd);cf=build();ci=build(True);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);rows=[]
 for inv,c,initial,desired in [(False,cf,src,dst),(True,ci,dst,src)]:
  sim=km.Simulator(batch_size=shots,seed=142137+inv);sim.use_same_registers_as(c)
  for k,v in initial.items():sim.write_across_shots(c.register_data[k],bits(v))
  sim.do(c)
  for k,v in desired.items():assert np.array_equal(sim.read_across_shots(c.register_data[k]),bits(v)),k
  public={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any();rows.append(dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like));sim.do(cf if inv else ci);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any()
 r=dict(status='PASS',shots=2*shots,rows=rows);Path('results/raw_three_jet_step_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()
