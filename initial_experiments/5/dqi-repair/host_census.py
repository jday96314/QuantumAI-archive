"""Generate exact host batches and relate sparse model decodes to host masks.

Instrumented sampling is intentionally kept here; official verifier is unchanged.
Includes equality checks against sample_shots for each batch's full syndromes,
masks and phases, then writes sparse error truth for the independent C++ census.
"""
import argparse,json
from pathlib import Path
import numpy as np,kickmix as km
from dqi.verify_submission_dqi import FieldTables,optimal_dqi_weights,ProblemSpecDQI,sample_shots,_initial_seed
p=argparse.ArgumentParser();p.add_argument('--seed',type=int);p.add_argument('--shots',type=int,default=102400);p.add_argument('--assess',action='store_true');a=p.parse_args()
spec=ProblemSpecDQI.benchmark_instance(num_shots=a.shots);tables=FieldTables(km.GF2Field(12));prob=optimal_dqi_weights(4095,35,12,2016)**2
base=_initial_seed('4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0','') if a.seed is None else a.seed
seeds=np.random.SeedSequence(base,pool_size=8).generate_state(8*(a.shots//256),dtype=np.uint64).reshape(-1,2,4)
events={c['shot']:c for c in map(json.loads,Path('results/exception_cases.jsonl').read_text().splitlines())} if a.assess else {}
rows=np.arange(1,71)[:,None];phase_errors=0;fixed_phase_errors=0;events_seen=0;weights_hist=np.zeros(36,int);first_hist=np.zeros(36,int)
f=None if a.assess else open('parameters/host_cases.bin','wb')
for batch,s in enumerate(seeds):
 seed=int.from_bytes(s[0],'little');rng=np.random.default_rng(seed);weights=rng.choice(36,size=256,p=prob);M=rng.integers(0,4095,size=(256,4095),endpoint=True,dtype=np.int64)
 synd=np.zeros((256,70),dtype=np.int64);ph=np.zeros(256,dtype=bool);records=np.zeros((256,141),dtype=np.uint16)
 for i,w in enumerate(weights):
  pos=rng.choice(4095,size=int(w),replace=False);val=rng.integers(1,4095,size=int(w),endpoint=True,dtype=np.int64)
  synd[i]=np.bitwise_xor.reduce(tables.antilog[(tables.log[val][None,:]+rows*pos[None,:])%4095],axis=1)
  ph[i]=sum(int(x).bit_count() for x in val&M[i,pos])%2
  records[i,:70]=synd[i];records[i,70]=w;records[i,71:71+w]=pos;records[i,106:106+w]=val
  c=events.get(batch*256+i)
  if c:
   assert c['syndromes']==synd[i].tolist();assert c['positions']==pos.tolist();assert c['values']==val.tolist()
   got=sum((v&int(M[i,j])).bit_count() for j,v in c['decoded'])%2
   phase_errors+=bool(got!=ph[i]);
   if "classical_fixed_decoded" in c:fixed_phase_errors+=bool(sum((v&int(M[i,j])).bit_count() for j,v in c["classical_fixed_decoded"])%2!=ph[i]);events_seen+=1;first_hist[c['first_stage']]+=1;weights_hist[w]+=1
 if not a.assess:
  # Every host-generated case must match our instrumented sampling exactly.
  ref=sample_shots(spec,tables,np.random.default_rng(seed),256,prob)
  assert np.array_equal(ref[0],synd) and np.array_equal(ref[1],M) and np.array_equal(ref[2],ph)
  records.tofile(f)
if f:f.close()
report=dict(status='PASS',shots=a.shots,base_seed=str(base),mode='assess' if a.assess else 'generate',host_exact_sampling_checked=not a.assess,exception_events=events_seen,model_phase_errors=phase_errors,classical_splice_phase_errors=fixed_phase_errors,first_gap_stage_hist=first_hist.tolist(),event_weight_hist=weights_hist.tolist())
Path('results/host_census_'+('assessment' if a.assess else 'generation')+'.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
