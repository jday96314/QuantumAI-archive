"""Reproduce canonical verifier failures for algorithm debugging, without scoring."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
import dqi.verify_submission_dqi as v
from stream_model import decode,pub
from published_dqi import isomorphism
args=argparse.ArgumentParser();args.add_argument('--tag',default='stream3');args.add_argument('--result');args.add_argument('--limit',type=int);a=args.parse_args()
HERE=Path(__file__).resolve().parent;c=km.Circuit.from_file(str(HERE/f'candidates/{a.tag}_full.kmb'))
summary=json.loads((HERE/'results'/(a.result or f'{a.tag}_verify.json')).read_text())
spec=v.ProblemSpecDQI.benchmark_instance(num_shots=102400);N=256
sq,mb,junk=v.validate_registers(c,spec);sq=km.array(sq);mb=km.array(mb)
tables=v.FieldTables(km.GF2Field(12));probs=v.optimal_dqi_weights(4095,35,12,2016)**2
base_seed=int(summary['circuit_sha256'],16);starts=range(0,102400,N)
seeds=np.random.SeedSequence(base_seed,pool_size=8).generate_state(8*len(starts),dtype=np.uint64).reshape(-1,2,4)
tasks=[(s,int.from_bytes(seeds[i,0],'little'),int.from_bytes(seeds[i,1],'little')) for i,s in enumerate(starts)]
def batch(task):
 start,rs,ss=task;s,m,ph=v.sample_shots(spec,tables,np.random.default_rng(rs),N,probs)
 sim=km.Simulator(batch_size=N,seed=ss);sim.use_same_registers_as(c)
 inputs=v.to_bits(s,12);sim.write_across_shots(sq,inputs);sim.write_across_shots(mb,v.to_bits(m,12));sim.do(c)
 ve=(sim.read_across_shots(sq)^inputs).any(0)|sim.read_across_shots(junk).any(0)
 pe=sim.read_phase_flipped_across_shots()^ph;ix=np.flatnonzero(ve|pe)
 return [(start+int(i),s[i],m[i],bool(ph[i]),bool(ve[i])) for i in ix]
bad=[x for group in v._run_parallel(batch,tasks,8) for x in group]
np.savez(HERE/f'results/{a.tag}_failures.npz',syndromes=np.array([x[1] for x in bad]),
         masks=np.array([x[2] for x in bad]),phases=np.array([x[3] for x in bad]))
f,back,tr=isomorphism(pub);counts={}
for index,ss,mm,phase,ve in (bad[:a.limit] if a.limit else bad):
 sig,why=decode([tr(int(x),f) for x in ss],offset=4 if 'offset4' in a.tag else 5);key=str(why);counts[key]=counts.get(key,0)+1
 print(json.dumps(dict(index=index,reason=why,value_error=ve)),flush=True)
print(json.dumps(dict(failures=len(bad),reasons=counts)),flush=True)
