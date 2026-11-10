"""Native circuit/model check on fixed boundary/random inputs and fresh seeds."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt.verify_submission_invsqrt import half_turns_to_units
from model import HERE,MAX_X,evaluate

def check(path,shots=65536,salt='phase-development-v1'):
    c=km.Circuit.from_file(str(path));rng=np.random.default_rng(7138801)
    boundary={1,MAX_X,958102}
    for k in range(22):
        for value in ((1<<k),3*(1<<k)):
            for d in range(-3,4):
                if 1<=value+d<=MAX_X:boundary.add(value+d)
    # Every normalization shift and every coefficient-prefix boundary.
    for k in range(22):
        for q in range(257):
            value=((1<<21)+(q<<13))>>k
            for d in (-1,0,1):
                if 1<=value+d<=MAX_X:boundary.add(value+d)
    base=np.array(sorted(boundary),dtype=np.int64)
    xs=np.resize(np.concatenate((base,np.arange(1,4097),rng.integers(1,MAX_X+1,max(0,shots-len(base)-4096)))) ,shots)
    expected=evaluate(xs);xreg=km.array(c.register_data['x']);public=set(xreg)
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    seed=int(hashlib.sha256((hashlib.sha256(path.read_bytes()).hexdigest()+salt).encode()).hexdigest(),16)
    batch=256;sim=km.Simulator(batch_size=batch,seed=seed,counting_operations=True);sim.use_same_registers_as(c)
    start=time.monotonic()
    for first in range(0,shots,batch):
        values=xs[first:first+batch];sim.clear_for_shot()
        bits=((values[None,:]>>np.arange(22)[:,None])&1).astype(bool)
        sim.write_across_shots(xreg,bits);sim.do(c)
        assert np.array_equal(sim.read_across_shots(xreg),bits),'input restoration'
        assert not sim.read_across_shots(scratch).any(),'dirty scratch'
        actual=np.array([half_turns_to_units(sim.read_shot_phase(j)) for j in range(batch)],dtype=np.uint64)
        bad=np.flatnonzero(actual!=expected[first:first+batch])
        assert not len(bad),(int(values[bad[0]]),int(actual[bad[0]]),int(expected[first+bad[0]]))
    stats=sim.operation_statistics()
    report=dict(status='PASS',scope='native simulation versus independent integer model; boundary/random subset',
        shots=shots,unique_inputs=int(len(np.unique(xs))),salt=salt,seconds=time.monotonic()-start,
        Q=c.num_qubits,mean_Toffolis=stats.num_toffoli_like/shots,
        rotations=stats.num_z_pow(min_clifford_hierarchy_level=4)/shots,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (HERE/'results'/f'{path.stem}_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate',type=Path);p.add_argument('--shots',type=int,default=65536)
    p.add_argument('--salt',default='phase-development-v1');a=p.parse_args();check(a.candidate,a.shots,a.salt)
