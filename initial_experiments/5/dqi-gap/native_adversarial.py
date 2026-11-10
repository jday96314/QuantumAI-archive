"""Native phase, restoration, scratch and repeated execution controls."""
import json
from pathlib import Path
import numpy as np,kickmix as km
from dqi.verify_submission_dqi import validate_registers,ProblemSpecDQI,to_bits
from classical_repair import divide,trim,product,plus
rng=np.random.default_rng(998231);cases=list(map(json.loads,Path('parameters/exception_cases_1m.jsonl').read_text().splitlines()));N=1024
# 768 exception cases, then 256 controls: zero, one error, weight 35 and
# random weights. Control syndromes must have a generic exact EEA path.
S=[];want=[];M=rng.integers(0,4096,(N,4095),dtype=np.int64);groups=[]
for i in range(N):
 if i<768:
  c=cases[i];s=c['syndromes'];decoded=c['decoded'];group='exception'
 else:
  kind=(i-768)%4;w=[0,1,35,int(rng.integers(2,36))][kind];group=['zero','one','weight35','random'][kind]
  while True:
   pos=rng.choice(4095,size=w,replace=False);val=rng.integers(1,4096,w);gf=km.GF2Field(12);s=[0]*70
   for j,v in zip(pos,val):
    for r in range(70):s[r]^=gf.mul(int(v),gf.pow(2,(r+1)*int(j)))
   a=[0]*70+[1];r=trim(s);fp=[];f=[1];good=True
   while r:
    q,rr=divide(a,r)
    if len(q)!=2:good=False;break
    fp,f=f,plus(fp,product(q,f));a,r=r,rr
    if len(f)-1>=w:break
   if good:break
  decoded=list(zip(pos.tolist(),val.tolist()))
 S.append(s);want.append(sum((int(v)&int(M[i,j])).bit_count() for j,v in decoded)%2);groups.append(group)
S=np.array(S,dtype=np.int64);want=np.array(want,dtype=bool);c=km.Circuit.from_file('../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb');syn,mask,anc=validate_registers(c,ProblemSpecDQI.benchmark_instance());sim=km.Simulator(batch_size=N,seed=791283);sim.use_same_registers_as(c);sim.write_across_shots(km.array(syn),to_bits(S,12));sim.write_across_shots(km.array(mask),to_bits(M,12));sim.do(c)
model_bad=np.flatnonzero(sim.read_phase_flipped_across_shots()^want).tolist();restored=bool(np.array_equal(sim.read_across_shots(km.array(syn)),to_bits(S,12)));clean=bool(not sim.read_across_shots(anc).any());sim.do(c)
result=dict(status='PASS' if not model_bad and restored and clean and not sim.read_phase_flipped_across_shots().any() else 'FAIL',shots=N,group_counts={g:groups.count(g) for g in set(groups)},model_phase_disagreements=model_bad,syndrome_restored=restored,scratch_clean=clean,twice_phase_clean=bool(not sim.read_phase_flipped_across_shots().any()),twice_syndrome_restored=bool(np.array_equal(sim.read_across_shots(km.array(syn)),to_bits(S,12))),twice_scratch_clean=bool(not sim.read_across_shots(anc).any()),scope='selected baseline simulator cross-check and adversarial controls, not acceptance of a replacement')
Path('results/native_adversarial.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
