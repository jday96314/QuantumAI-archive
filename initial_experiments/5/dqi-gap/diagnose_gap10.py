"""Matched inputs for the full forced-gap diagnostic and selected baseline.

These shots are NOT drawn from the acceptance distribution. No R or DQI cost
is inferred from their failure rate. Includes out-of-domain controls explicitly
to reject the forced branch as a standalone replacement.
"""
import json
from pathlib import Path
import numpy as np,kickmix as km
from dqi.verify_submission_dqi import validate_registers,ProblemSpecDQI,to_bits
from classical_repair import mul
rng=np.random.default_rng(2026102405);pool=[c for c in map(json.loads,Path('parameters/exception_cases_1m.jsonl').read_text().splitlines()) if c['first_stage']==10 and c['gap_count']==1]
assert pool
N=512;S=[];truth=[];M=rng.integers(0,4096,(N,4095),dtype=np.int64);group=[]
for i in range(N):
 if i<448:
  c=pool[i%len(pool)];s=c['syndromes'];pos=c['positions'];val=c['values'];g='gap10'
 else:
  kind=(i-448)//16;w=[0,1,35,12][kind];g=['zero','one','weight35','weight12'][kind];pos=rng.choice(4095,w,replace=False).tolist();val=rng.integers(1,4096,w).tolist();s=[0]*70
  gf=km.GF2Field(12)
  for j,v in zip(pos,val):
   for r in range(70):s[r]^=gf.mul(v,gf.pow(2,(r+1)*j))
 S.append(s);truth.append(sum((v&int(M[i,j])).bit_count() for j,v in zip(pos,val))%2);group.append(g)
S=np.array(S,dtype=np.int64);truth=np.array(truth,dtype=bool);np.savez('parameters/gap10_matched_shots.npz',syndromes=S,mbuc=M,phases=truth,groups=np.array(group));reports=[]
for name,path in [('selected','../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb'),('gap10_forced','candidates/gap10_diagnostic.kmb')]:
 c=km.Circuit.from_file(path);syn,mask,anc=validate_registers(c,ProblemSpecDQI.benchmark_instance());sim=km.Simulator(batch_size=N,seed=813829,counting_operations=True);sim.use_same_registers_as(c);sim.write_across_shots(km.array(syn),to_bits(S,12));sim.write_across_shots(km.array(mask),to_bits(M,12));sim.do(c)
 vbad=(sim.read_across_shots(km.array(syn))^to_bits(S,12)).any(0)|sim.read_across_shots(anc).any(0);pbad=sim.read_phase_flipped_across_shots()^truth;bygroup={}
 for g in sorted(set(group)):
  take=np.array(group)==g;bygroup[g]=dict(shots=int(take.sum()),value_errors=int((vbad&take).sum()),phase_errors=int((pbad&take&~vbad).sum()))
 T=sim.operation_statistics().num_toffoli_like/N;sim.do(c);take=np.array(group)=='gap10';twicebad=(sim.read_across_shots(km.array(syn))^to_bits(S,12)).any(0)|sim.read_across_shots(anc).any(0)|sim.read_phase_flipped_across_shots();reports.append(dict(name=name,Q=c.num_qubits,executed_mean_Toffolis=T,static_Toffolis=c.operation_statistics().num_toffoli_like,groups=bygroup,twice_gap10_errors=int((twicebad&take).sum())))
result=dict(status='DIAGNOSTIC_ONLY',seed=2026102405,unique_gap10_cases=len(pool),scope=__doc__,circuits=reports)
Path('results/gap10_native_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
