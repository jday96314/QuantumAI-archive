"""Native quotient test, including both denominator axes and zero."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from extension_inverse import Adapter,setup
from schur_dynamic_division_sequential import create
m,k,fw,bk,tr=setup();b=FastBuilder(km.CircuitBuilder());num=b.create_quantum_register(24,'num');den=b.create_quantum_register(24,'den');out=b.create_quantum_register(24,'out');flags=b.create_quantum_register(2,'flags');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));deltas=[m.gmul(tr(2,fw),tr(2,fw))^tr(2,fw),m.gmul(tr(6,fw),tr(6,fw))^tr(6,fw)];create(b,E,[q.id for q in num],[q.id for q in den],[q.id for q in out],k,list(zip([q.id for q in flags],deltas)));ab.release_recycle('test',sc);c=b.finish_circuit();c=km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)));shots=16384;rng=np.random.default_rng(142099);aa=rng.integers(0,1<<24,shots);dd=rng.integers(0,1<<24,shots);dd[:4096]=np.arange(4096);dd[4096:8192]=np.arange(4096)<<12;ff=rng.integers(0,4,shots)
def div(x,y,k):
 a,b=x&4095,x>>12;d,e=y&4095,y>>12;nn=m.gmul(d,d)^m.gmul(d,e)^m.gmul(k,m.gmul(e,e));ni=m.ginv(nn);u=m.gmul(m.gmul(a,d^e)^m.gmul(k,m.gmul(b,e)),ni);t=m.gmul(m.gmul(a,e)^m.gmul(b,d),ni);return u|(t<<12)
want=np.array([div(int(x),int(y),k^(deltas[0] if f&1 else 0)^(deltas[1] if f&2 else 0)) for x,y,f in zip(aa,dd,ff)]);bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142099);sim.use_same_registers_as(c)
for key,v in [('num',aa),('den',dd),('flags',ff)]:sim.write_across_shots(c.register_data[key],bits(v,len(c.register_data[key])))
sim.do(c);errors=[]
for key,v in [('num',aa),('den',dd),('flags',ff),('out',want)]:
 if not np.array_equal(sim.read_across_shots(c.register_data[key]),bits(v,len(c.register_data[key]))):errors.append(key)
public={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);r=dict(status='PASS' if not errors and not sim.read_across_shots(junk).any() and not sim.read_phase_flipped_across_shots().any() else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=shots,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));print(json.dumps(r),flush=True);Path('results/sequential_division_native.json').write_text(json.dumps(r,indent=2)+'\n');assert r['status']=='PASS'
