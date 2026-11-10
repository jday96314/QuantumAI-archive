"""Check that intentionally unrepaired residual paths still restore all data."""
import json
from pathlib import Path
import numpy as np,kickmix as km
from dqi.verify_submission_dqi import validate_registers,ProblemSpecDQI,to_bits
from classical_repair import divide,trim,product,plus

cases=list(map(json.loads,Path('parameters/residual_inputs.jsonl').read_text().splitlines()))
assert [c['weight'] for c in cases]==[29,32,33,30,30]
per=128;N=per*len(cases);rng=np.random.default_rng(2026102507);M=rng.integers(0,4096,(N,4095),dtype=np.int64)
S=np.asarray([c['syndromes'] for c in cases for _ in range(per)],dtype=np.int64)
want=np.asarray([sum((v&int(M[i,j])).bit_count() for j,v in zip(cases[i//per]['positions'],cases[i//per]['values']))%2 for i in range(N)],dtype=bool)
manifest=json.loads(Path('parameters/best_candidate_manifest.json').read_text());c=km.Circuit.from_file(manifest['circuit'])
syn,mask,anc=validate_registers(c,ProblemSpecDQI.benchmark_instance());sim=km.Simulator(batch_size=N,seed=710521);sim.use_same_registers_as(c)
sim.write_across_shots(km.array(syn),to_bits(S,12));sim.write_across_shots(km.array(mask),to_bits(M,12));sim.do(c)
phase=sim.read_phase_flipped_across_shots()^want;value=(sim.read_across_shots(km.array(syn))!=to_bits(S,12)).any(0)|sim.read_across_shots(anc).any(0)
rows=[]
for i,cc in enumerate(cases):
    P=[0]*70+[1];R=trim(cc['syndromes']);bp=[];b=[1];pattern=[]
    while R:
        q,r=divide(P,R);bp,b=b,plus(bp,product(q,b));P,R=R,r;pattern.append(len(q)-1)
        if len(b)-1>=cc['weight']:break
    rows.append(dict(seed=cc['seed'],shot=cc['shot'],weight=cc['weight'],quotient_degrees=pattern,mask_cases=per,truth_phase_errors=int(phase[i*per:(i+1)*per].sum()),value_errors=int(value[i*per:(i+1)*per].sum())))
sim.do(c);twice=dict(phase_errors=int(sim.read_phase_flipped_across_shots().sum()),value_errors=int(((sim.read_across_shots(km.array(syn))!=to_bits(S,12)).any(0)|sim.read_across_shots(anc).any(0)).sum()))
report=dict(status='PASS' if not value.any() and not any(twice.values()) else 'FAIL',scope='known classical residuals: phase mistakes are retained, syndrome/scratch preservation and twice-execution phase restoration required',shots=N,cases=rows,twice=twice)
Path('results/best_residual_preservation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
raise SystemExit(0 if report['status']=='PASS' else 1)
