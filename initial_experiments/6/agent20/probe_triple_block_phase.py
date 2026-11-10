"""Locate the first nonzero phase in the complete triple-block composition."""
from pathlib import Path
import inspect,json,subprocess
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
from compact_blocked_cfraction import Arithmetic
import triple_gap_block as block
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent
m=load();m.inv_word=bind(m);b=m.Builder();N=8;W=[b.qubits(12) for _ in range(N)];E=m.Em(b,*b.qubits(3),m.Alloc(b));marks=[]
source=inspect.getsource(block.block)
replacements={
 'copy();head.apply(m,E,H)':'copy();head.apply(m,E,H);mark("head")',
 '   census()\n   with frame.parameters':'   census();mark("census")\n   with frame.parameters',
 '    with ExitStack() as stack:':'    mark("frame_enter")\n    with ExitStack() as stack:',
 'pr);E.mac(pr)':'pr);E.mac(pr)',
 'pr=m.polyprod_rt(A[1:],raw,targets);E.mac(pr)':'pr=m.polyprod_rt(A[1:],raw,targets);E.mac(pr);mark("P_prepare")',
 '     E.free(pr)':'     mark("tail_done");E.free(pr);mark("P_free")',
 '   census()\n  head.apply':'   mark("frame_out");census();mark("census_out")\n  head.apply',
 '  head.apply(m,E,H,True);copy()':'  head.apply(m,E,H,True);mark("head_inverse");copy();mark("copy_clear")',
 '  if not inverse:head.apply(m,E,H)':'  if not inverse:head.apply(m,E,H);mark("head_final")',
}
for old,new in replacements.items():
 if old==new:continue
 assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
(HERE/'parameters/triple_gap_block_phase_probe_source.py').write_text(source)
ns=dict(block.__dict__,mark=lambda name:marks.append((len(b._items),name)));exec(source,ns);ns['block'](m,E,Arithmetic(m,E),W)
rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tape=HERE/'parameters/triple_gap_block_phase_probe.npy';np.save(tape,rows);kb=km.CircuitBuilder()
for j in range(N):kb.create_quantum_register(12,f'W{j}')
pre=HERE/'parameters/triple_gap_block_phase_probe.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=HERE/'parameters/empty';rng=np.random.default_rng(2026102608);shots=512;values=rng.integers(0,4096,(shots,N));values[:3]=[[0]*N,[1]*N,[1,0,*([0]*(N-2))]];values[3,:4]=[0,1,2,3];values[4,:4]=[1,0,2,3];bits=lambda x:((x[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
for stop,name in marks:
 out=HERE/'candidates'/f'triple_phase_probe_{name}.kmx';subprocess.run([str(ref.probe_translator()),str(tape),str(stop),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(stop)],check=True,capture_output=True);c=km.Circuit.from_file(str(out),format='kmx');sim=km.Simulator(batch_size=shots,seed=20261026);sim.use_same_registers_as(c)
 for j in range(N):sim.write_across_shots(c.register_data[f'W{j}'],bits(values[:,j]))
 sim.do(c);phase=sim.read_phase_flipped_across_shots();r=dict(stage=name,phase_errors=int(phase.sum()),indices=np.flatnonzero(phase).tolist(),T=c.operation_statistics().num_toffoli_like,Q=c.num_qubits);reports.append(r);print(json.dumps({k:v for k,v in r.items() if k!='indices'}),flush=True)
(HERE/'results/triple_block_phase_probe.json').write_text(json.dumps(reports,indent=2)+'\n')
