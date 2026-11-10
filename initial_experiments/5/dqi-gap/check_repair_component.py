"""Rebuild and independently validate reversible local Padé exception repair."""
import json,time
from pathlib import Path
import numpy as np,kickmix as km
import lookahead_six_quantum as q
from import_published import load
from low_workspace_inverse import bind
from published_dqi import isomorphism
from classical_repair import repair,mul
# All outputs from the reused build routine are redirected into this experiment.
q.ROOT=Path.cwd()
# The translator remains read-only; a local symlink keeps q.build self-contained.
link=Path('translate_prefix')
if not link.exists():link.symlink_to('../../4/codex/translate_prefix')
rng=np.random.default_rng(819213);shots=4096;m=load();fw,bw,tr=isomorphism(m)
rows=[];wants=[];singular=0
while len(rows)<shots:
 F=rng.integers(0,4096,6).tolist()
 if len(rows)%2:F[0]=0
 try:U,V,pi=repair(F)
 except ValueError:singular+=1;continue
 # Circuit's physical order is v1,v0,u0,u1,pinv,v2.
 rows.append([tr(x,fw) for x in F]);wants.append([tr(x,fw) for x in [V[1],V[0],U[0],U[1],pi,V[2]]])
rows=np.array(rows).T;wants=np.array(wants).T
bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
reports=[]
for which in ['forward','inverse','both']:
 c=q.build(which);sim=km.Simulator(batch_size=shots,seed=991289);sim.use_same_registers_as(c)
 start=wants if which=='inverse' else rows;want=wants if which=='forward' else rows
 for name,v in zip(['e','d','a','b','c','f'],start):sim.write_across_shots(c.register_data[name],bits(v))
 sim.do(c)
 for name,v in zip(['e','d','a','b','c','f'],want):assert np.array_equal(sim.read_across_shots(c.register_data[name]),bits(v)),(which,name)
 public={qq for reg in c.register_data.values() for qq in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
 assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
 reports.append(dict(direction=which,shots=shots,Q=c.num_qubits,static_Toffolis=c.operation_statistics().num_toffoli_like))
result=dict(status='PASS',scope='reversible six-moment correction matched independent native-field polynomial algebra; 2048 forced-zero first pivots in each direction; NOT complete decoder',circuits=reports,singular_rejected=singular)
Path('results/repair_component.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
