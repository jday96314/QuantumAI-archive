"""Common-case regression for saved failures, without resampling for a score."""
import json,sys
from pathlib import Path
import numpy as np
import kickmix as km
files=sorted(Path('results/counterfree_failures').glob('*.npz'));zs=[dict(np.load(p)) for p in files];ss=np.concatenate([z['syndromes'] for z in zs]);ee=np.concatenate([z['mbuc'] for z in zs]);pp=np.concatenate([z['phases'] for z in zs]);ids=np.concatenate([z['shot_indices'] for z in zs]);bits=lambda v:((v[:,None,:]>>np.arange(12)[None,:,None])&1).reshape(-1,len(v)).astype(bool)
# Explicit register writes avoid assumptions about serialized register ordering.
rows=[]
for path in sys.argv[1:]:
 c=km.Circuit.from_file(path,format='kmb');source=km.array([q for i in range(70) for q in c.register_data[f'S{i}']]);cl=km.array([q for i in range(4095) for q in c.register_data[f'E{i}']]);public=set(source);junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public])
 for seed in [142103,142104]:
  sim=km.Simulator(batch_size=len(ss),seed=seed);sim.use_same_registers_as(c)
  for i in range(70):sim.write_across_shots(c.register_data[f'S{i}'],((ss[:,i][None,:]>>np.arange(12)[:,None])&1).astype(bool))
  for i in range(4095):sim.write_across_shots(c.register_data[f'E{i}'],((ee[:,i][None,:]>>np.arange(12)[:,None])&1).astype(bool))
  initial=sim.read_across_shots(source);sim.do(c);bad=(sim.read_across_shots(source)!=initial).any(0)|sim.read_across_shots(junk).any(0);phase=(sim.read_phase_flipped_across_shots()!=pp)&~bad;r=dict(path=path,seed=seed,shots=len(ss),value_errors=ids[bad].tolist(),phase_errors=ids[phase].tolist());rows.append(r);print(json.dumps(r),flush=True)
Path('results/counterfree_failures_diagnosis.json').write_text(json.dumps(rows,indent=2)+'\n')
