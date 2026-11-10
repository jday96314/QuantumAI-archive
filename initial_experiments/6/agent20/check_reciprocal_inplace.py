"""Exhaust all field inputs, source preservation of the inverse sandwich, phase."""
from pathlib import Path
import json,subprocess,hashlib
import numpy as np
import kickmix as km
from import_published import load
import reciprocal_inplace as reciprocal
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent


def run(controlled=False):
    m=load();b=m.Builder();src=b.qubits(12);control=b.qubits(1) if controlled else []
    E=m.Em(b,*b.qubits(3),m.Alloc(b));reciprocal.apply(m,E,src,control[0] if controlled else None)
    tag='reciprocal_inplace'+('_controlled' if controlled else '')
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder();kb.create_quantum_register(12,'x')
    if controlled:kb.create_quantum_register(1,'c')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb')
    values=np.tile(np.arange(4096),2 if controlled else 1);flags=np.repeat([False,True],4096) if controlled else None
    want=np.array([m.ginv(int(v)) for v in values]);want=np.where(flags,want,values) if controlled else want
    bits=lambda x:((np.asarray(x)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    public=12+int(controlled)
    for seed in (2026102508,2026102509):
        sim=km.Simulator(batch_size=len(values),seed=seed);sim.use_same_registers_as(c)
        sim.write_across_shots(c.register_data['x'],bits(values))
        if controlled:sim.write_across_shots(c.register_data['c'],flags[None,:])
        sim.do(c);assert np.array_equal(sim.read_across_shots(c.register_data['x']),bits(want))
        if controlled:assert np.array_equal(sim.read_across_shots(c.register_data['c']),flags[None,:])
        assert not sim.read_phase_flipped_across_shots().any(),'phase'
        assert not sim.read_across_shots(km.array([km.q(i) for i in range(public,c.num_qubits)])).any(),'scratch'
        sim.do(c);assert np.array_equal(sim.read_across_shots(c.register_data['x']),bits(values))
        assert not sim.read_phase_flipped_across_shots().any(),'twice phase'
        assert not sim.read_across_shots(km.array([km.q(i) for i in range(public,c.num_qubits)])).any(),'twice scratch'
    row=dict(status='PASS_EXHAUSTIVE_COMPONENT',controlled=controlled,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,
             shots=2*len(values),sha256=hashlib.sha256(out.with_suffix('.kmb').read_bytes()).hexdigest(),
             scope='all4096 field values including zero; both controls if present; two measurement seeds; involution and clean scratch/phase')
    print(json.dumps(row),flush=True);return row

if __name__=='__main__':
    rows=[run(),run(True)];(HERE/'results/reciprocal_inplace_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
