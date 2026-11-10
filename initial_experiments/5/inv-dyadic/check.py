"""Native whole-domain random and prefix-boundary arithmetic/phase checks."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from certify import evaluate,reference,MAX_X,HERE

p=argparse.ArgumentParser();p.add_argument('circuit',type=Path);p.add_argument('model',type=Path)
p.add_argument('--random',type=int,default=32768);p.add_argument('--repeats',type=int,default=4);a=p.parse_args()
c=km.Circuit.from_file(str(a.circuit));m=json.loads(a.model.read_text());wc=m['widths'][0]
rng=np.random.default_rng(20261024);cases={1,2,3,MAX_X,2606235,958599}
for pp,lo,d,f,degree,*rest in m['leaves']:
 for k in range(22):
    if k%2!=pp:continue
    for n in [2**21+lo,2**21+lo+2**(21-d)-1]:
        q=n>>k
        cases.update(x for x in [q-1,q,q+1] if 1<=x<=MAX_X)
xs=np.concatenate([np.array(sorted(cases)),rng.integers(1,MAX_X+1,a.random)])
xs=np.resize(xs,((len(xs)+255)//256)*256);expected=evaluate(m,xs).astype(np.uint64)<<np.uint64(64-wc-10)
xreg=c.register_data['x'];scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in xreg])
target,tol=reference();t0=time.time();peak_error=0.;total=0;mean_t=0
for repeat in range(a.repeats):
 sim=km.Simulator(batch_size=256,seed=7171+repeat,counting_operations=True);sim.use_same_registers_as(c)
 for start in range(0,len(xs),256):
    xx=xs[start:start+256];sim.clear_for_shot()
    bits=((xx[None,:]>>np.arange(22)[:,None])&1).astype(bool)
    sim.write_across_shots(xreg,bits);sim.do(c)
    assert np.array_equal(sim.read_across_shots(xreg),bits),'input not restored'
    assert not sim.read_across_shots(scratch).any(),'dirty scratch'
    phase=[]
    for i in range(256):
        n,d=sim.read_shot_phase(i).as_integer_ratio();phase.append((((n<<64)//d+1)>>1)&((1<<64)-1))
    phase=np.asarray(phase,dtype=np.uint64)
    assert np.array_equal(phase,expected[start:start+256]),f'phase differs from integer model near {xx[0]}'
    used=np.abs((phase-target[xx-1]).view(np.int64)).astype(np.uint64)
    assert (used<tol[xx-1]).all(),'phase tolerance failed'
    peak_error=max(peak_error,float((used/tol[xx-1]).max()));total+=256
 mean_t+=sim.operation_statistics().num_toffoli_like
report=dict(status='PASS',scope='random and prefix boundaries; not full host verification',circuit=str(a.circuit),
 model=str(a.model),verified_shots=total,repeats=a.repeats,worst_error_over_tolerance=peak_error,
 mean_toffolis=mean_t/total,Q=c.num_qubits,seconds=time.time()-t0)
(HERE/'results'/f'{a.circuit.stem}_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
