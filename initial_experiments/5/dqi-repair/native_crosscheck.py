"""Compare independent field-level compact predictions against native phase."""
import json,time
from pathlib import Path
import numpy as np,kickmix as km
from dqi.verify_submission_dqi import verify,ProblemSpecDQI
p=Path('results/exception_cases.jsonl'); cases=[json.loads(s) for s in p.read_text().splitlines()]
rng=np.random.default_rng(2718281828);N=256
S=np.array([cases[i%len(cases)]['syndromes'] for i in range(N)],dtype=np.int64)
M=rng.integers(0,4096,(N,4095),dtype=np.int64)
expected=np.zeros(N,dtype=bool);truth=np.zeros(N,dtype=bool)
for i in range(N):
 case=cases[i%len(cases)]
 expected[i]=sum((v&int(M[i,j])).bit_count() for j,v in case['decoded'])%2
 truth[i]=sum((v&int(M[i,j])).bit_count() for j,v in zip(case['positions'],case['values']))%2
np.savez('parameters/native_crosscheck_shots.npz',syndromes=S,mbuc=M,predicted=expected,truth=truth)
c=km.Circuit.from_file('../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb')
# Direct simulation allows census/native disagreement rates above host acceptance.
from dqi.verify_submission_dqi import validate_registers,to_bits
syn,mask,anc=validate_registers(c,ProblemSpecDQI.benchmark_instance())
sim=km.Simulator(batch_size=N,seed=992912,counting_operations=True);sim.use_same_registers_as(c)
sim.write_across_shots(km.array(syn),to_bits(S,12));sim.write_across_shots(km.array(mask),to_bits(M,12));sim.do(c)
phase=sim.read_phase_flipped_across_shots();bad=np.flatnonzero(phase^expected).tolist()
result=dict(cases=N,unique_exceptions=len(cases),native_model_phase_disagreements=len(bad),disagreement_indices=bad,native_truth_phase_errors=int((phase^truth).sum()),syndrome_restored=bool(np.array_equal(sim.read_across_shots(km.array(syn)),to_bits(S,12))),scratch_clean=bool(not sim.read_across_shots(anc).any()))
Path('results/native_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
