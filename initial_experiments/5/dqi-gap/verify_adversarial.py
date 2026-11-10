"""Actual native phase contract on independently masked exceptions and controls."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,kickmix as km
from dqi.verify_submission_dqi import validate_registers,ProblemSpecDQI,to_bits
from serial_reference import decode

p=argparse.ArgumentParser();p.add_argument('circuit');p.add_argument('--shots',type=int,default=1024);p.add_argument('--m',type=int,default=4095);p.add_argument('--seed',type=int,default=2026102501);p.add_argument('--output',default='results/adversarial.json');a=p.parse_args()
rng=np.random.default_rng(a.seed);cases=list(map(json.loads,Path('parameters/exception_cases_1m.jsonl').read_text().splitlines()));N=a.shots
preferred=[]
for stage in range(34):
    preferred += [c for c in cases if c['first_stage']==stage][:1]
preferred += [c for c in cases if decode(c['syndromes'])[-1]['terminal'] or c['gap_count']>1]
seen=set();ordered=[]
for cc in preferred+cases:
    if cc['shot'] not in seen:ordered.append(cc);seen.add(cc['shot'])
cases=ordered
S=[];want=[];groups=[];M=rng.integers(0,4096,(N,a.m),dtype=np.int64);gf=km.GF2Field(12)
for i in range(N):
    if i<3*N//4:
        cc=cases[i%len(cases)];s=cc['syndromes'];pos=cc['positions'];val=cc['values'];group=f'gap_stage_{cc["first_stage"]}'
        meta=decode(s)[-1]
        if meta['terminal']:group+=':immediate_completion'
        if cc['gap_count']>1:group+=':multiple'
    else:
        kind=(i-3*N//4)%4;w=[0,1,35,int(rng.integers(2,36))][kind];group=['zero','one','weight35','random'][kind]
        pos=rng.choice(4095,size=w,replace=False).tolist();val=rng.integers(1,4096,w).tolist();s=[0]*70
        for j,v in zip(pos,val):
            for r in range(70):s[r]^=gf.mul(v,gf.pow(2,(r+1)*j))
    S.append(s);want.append(sum((v&int(M[i,j])).bit_count() for j,v in zip(pos,val) if j<a.m)%2);groups.append(group)
S=np.asarray(S,dtype=np.int64);want=np.asarray(want,dtype=bool);c=km.Circuit.from_file(a.circuit);syn,mask,anc=validate_registers(c,ProblemSpecDQI.benchmark_instance(m=a.m))
sim=km.Simulator(batch_size=N,seed=a.seed^0x798572,counting_operations=True);sim.use_same_registers_as(c)
sim.write_across_shots(km.array(syn),to_bits(S,12));sim.write_across_shots(km.array(mask),to_bits(M,12));sim.do(c)
phase=sim.read_phase_flipped_across_shots()^want;value=(sim.read_across_shots(km.array(syn))!=to_bits(S,12)).any(0)|sim.read_across_shots(anc).any(0)
stats=sim.operation_statistics();rows={}
for group in sorted(set(groups)):
    use=np.array([x==group for x in groups]);rows[group]=dict(shots=int(use.sum()),phase_errors=int(phase[use].sum()),value_errors=int(value[use].sum()))
sim.do(c);twice=dict(phase_errors=int(sim.read_phase_flipped_across_shots().sum()),value_errors=int(((sim.read_across_shots(km.array(syn))!=to_bits(S,12)).any(0)|sim.read_across_shots(anc).any(0)).sum()))
report=dict(status='PASS' if not phase.any() and not value.any() and not any(twice.values()) else 'FAIL',scope='adversarial diagnostic; no DQI score inferred',shots=N,seed=a.seed,circuit=a.circuit,sha256=hashlib.sha256(Path(a.circuit).read_bytes()).hexdigest(),Q=c.num_qubits,mean_executed_Toffolis=int(stats.num_toffoli_like)/N,static_Toffolis=c.operation_statistics().num_toffoli_like,phase_errors=int(phase.sum()),value_errors=int(value.sum()),twice=twice,groups=rows)
Path(a.output).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
raise SystemExit(0 if report['status']=='PASS' else 1)
