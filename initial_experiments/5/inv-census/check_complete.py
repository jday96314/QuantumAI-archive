"""Frozen complete-circuit model matches; no scoring/verification replacement."""
import hashlib,json,math
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt.verify_submission_invsqrt import validate_registers,ProblemSpecInvSqrt,to_bits,from_bits,half_turns_to_units
from prototype import evaluate,N
H=Path(__file__).resolve().parent

def check(mode):
    path=H/f'candidates/{mode}_best.kmb';c=km.Circuit.from_file(str(path));digest=hashlib.sha256(path.read_bytes()).hexdigest()
    x,scratch=validate_registers(c,ProblemSpecInvSqrt.for_n_p(10))
    rng=np.random.default_rng(261024)
    vals=[1,2,3,4,958102,2965859,N]
    for k in range(22):
        for a in range(257):
            v=((256+a)*(1<<k))//256
            vals.extend(z for z in (v-1,v,v+1) if 1<=z<=N)
    xs=np.unique(np.r_[vals,rng.integers(1,N+1,16384)]).astype(np.int64)
    original=xs.copy();total=0;mean=0
    for salt in ('census-complete-independent-a','census-complete-independent-b'):
        count=math.ceil(len(xs)/256)*256;xxs=np.resize(xs,count)
        sim=km.Simulator(batch_size=256,counting_operations=True,seed=int(hashlib.sha256(salt.encode()).hexdigest(),16));sim.use_same_registers_as(c)
        for first in range(0,count,256):
            xx=xxs[first:first+256];sim.clear_for_shot();sim.write_across_shots(x,to_bits(xx.astype(np.uint64),22));sim.do(c)
            assert np.array_equal(from_bits(sim.read_across_shots(x)),xx)
            assert not sim.read_across_shots(scratch).any()
            k,m,s,p=evaluate(xx)
            expected=(p<<s).astype(np.uint64)<<np.uint64(27)
            actual=np.array([half_turns_to_units(sim.read_shot_phase(j)) for j in range(256)],dtype=np.uint64)
            assert np.array_equal(actual,expected),(mode,salt,int(xx[np.argmax(actual!=expected)]))
        total+=count;mean+=sim.operation_statistics().num_toffoli_like
    return dict(status='EXACT_NATIVE_MODEL_MATCH',scope='development only, not full native acceptance',
        sha256=digest,unique_inputs=len(original),shots=total,seeds=2,Q=c.num_qubits,mean_Toffolis=mean/total)

if __name__=='__main__':
    results={}
    for mode in ('borrowed','circular','walsh'):
        results[mode]=check(mode);print(json.dumps(results[mode]),flush=True)
    (H/'results/complete_model_matches.json').write_text(json.dumps(results,indent=2)+'\n')
