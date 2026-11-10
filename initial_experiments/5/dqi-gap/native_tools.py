"""Local tape translation and simulator checks; the host verifier is untouched."""
import json,subprocess,os
from pathlib import Path
import numpy as np,kickmix as km


def finish(b,registers,stem):
    stem += os.environ.get('DQI_GAP_TAG','')
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tape=Path('parameters')/(stem+'.npy');np.save(tape,rows)
    kb=km.CircuitBuilder()
    for name,word in registers:kb.create_quantum_register(len(word),name)
    pre=Path('parameters')/(stem+'.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=Path('parameters/empty');empty.write_text('');out=Path('candidates')/(stem+'.kmx')
    subprocess.run(['./translate_prefix',str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def check(c,initial,want,shots,seed=81276):
    bits=lambda v,n:((np.asarray(v)[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
    for name,values in initial.items():sim.write_across_shots(c.register_data[name],bits(values,len(c.register_data[name])))
    sim.do(c);errors={}
    for name,values in want.items():
        bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(values,len(c.register_data[name]))).any(0))
        if len(bad):errors[name]=bad[:12].tolist()
    public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots())
    return dict(status='PASS' if not errors and not len(junk) and not len(phase) else 'FAIL',errors=errors,junk=junk[:12].tolist(),phase=phase[:12].tolist(),Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)
