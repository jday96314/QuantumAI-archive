"""Compare converter charts on common failures captured before this fix."""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from derivative_inverse_shifted import build as oldbuild
from derivative_inverse_lean import build as newbuild
z=dict(np.load('parameters/counterfree_failure_polynomials.npz'));A=z['A'];B=z['B'];ids=z['shot_indices'];c=np.array([next((i for i,v in enumerate(b.tolist()+[1]) if v),35)+1 for b in B]);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);out=[]
for name,build in [('old',oldbuild),('new',newbuild)]:
 circ=build(35);sim=km.Simulator(batch_size=len(A),seed=142106);sim.use_same_registers_as(circ)
 for i in range(35):
  for key,values in [('A',A),('B',B)]:sim.write_across_shots(circ.register_data[f'{key}{i}'],bits(values[:,i]))
 if name=='old':sim.write_across_shots(circ.register_data['padding'],bits(c,len(circ.register_data['padding'])))
 sim.do(circ);aa=np.stack([sim.read_across_shots(circ.register_data[f'A{i}']) for i in range(35)]);pub={q for reg in circ.register_data.values() for q in reg};junk=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in pub]);r=dict(name=name,Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like,scratch=ids[sim.read_across_shots(junk).any(0)].tolist(),phase=ids[sim.read_phase_flipped_across_shots()].tolist());print(json.dumps(r),flush=True);out.append((r,aa))
bad=ids[(out[0][1]!=out[1][1]).any((0,1))].tolist();r=dict(status='PASS' if not bad and out[0][0]['scratch']==out[1][0]['scratch'] and out[0][0]['phase']==out[1][0]['phase'] else 'FAIL',shots=len(A),different_outputs=bad,rows=[x[0] for x in out]);Path('results/counterfree_fixed_regressions.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);assert r['status']=='PASS'
