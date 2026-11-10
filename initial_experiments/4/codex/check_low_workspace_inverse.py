"""Exhaustive field and random-target check of the reduced inversion lifetime."""
import json,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
ROOT=Path(__file__).resolve().parent
m=load();b=m.Builder();source=b.qubits(12);target=b.qubits(12);E=m.Em(b,*b.qubits(3),m.Alloc(b))
with bind(m)(E,source) as iv:
 for x,y in zip(iv,target):E.cx(x,y)
rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag='low_workspace_inverse';tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows)
kb=km.CircuitBuilder();kb.create_quantum_register(12,'source');kb.create_quantum_register(12,'target');pre=ROOT/f'parameters/{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{tag}.empty';empty.write_text('');full=ROOT/f'candidates/{tag}.kmx'
subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(full),str(b.next_out),str(empty),str(len(rows))],check=True)
c=km.Circuit.from_file(str(full),format='kmx');rng=np.random.default_rng(91025);values=np.arange(4096);incoming=rng.integers(0,4096,4096);want=incoming^np.array([m.ginv(int(v)) for v in values]);bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
s=km.Simulator(batch_size=4096,seed=18);s.use_same_registers_as(c);s.write_across_shots(c.register_data['source'],bits(values));s.write_across_shots(c.register_data['target'],bits(incoming));s.do(c)
assert np.array_equal(s.read_across_shots(c.register_data['source']),bits(values));assert np.array_equal(s.read_across_shots(c.register_data['target']),bits(want));assert not s.read_phase_flipped_across_shots().any()
public=set(c.register_data['source'])|set(c.register_data['target']);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not s.read_across_shots(scratch).any()
r=dict(status='PASS',scope='all4096 field values including zero, random incoming XOR targets, source/scratch restoration and phase',shots=4096,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)
Path('results/low_workspace_inverse_check.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
