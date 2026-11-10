"""Native value, phase, and scratch tests for a quantum shift selector."""
import json,subprocess
from pathlib import Path
import kickmix as km
import numpy as np
from import_published import load
from cube_arithmetic import CubeArithmetic
import check_blocked_cfraction as reference
from variable_shift_inverse import apply

HERE=Path(__file__).resolve().parent


def build(N):
    m=load();b=m.Builder();R=[b.qubits(12) for _ in range(N)];flag=b.qubits(1)[0]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=CubeArithmetic(m,E);apply(m,E,math,R,flag)
    tag=f'variable_shift_n{N}';tape=HERE/'parameters'/f'{tag}.npy'
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
    kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'R{i}')
    kb.create_quantum_register(1,'extra');pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');text=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(reference.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(text),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(text),format='kmx');c.to_file(str(text.with_suffix('.kmb')),format='kmb');return c


def check(N,shots=512):
    m=load();rng=np.random.default_rng(2026102800+N);values=rng.integers(0,4096,(shots,N)).tolist();flags=np.arange(shots)%2
    values[:2]=[[0]*N,[1]*N]
    expected=[reference.divide(m,row,[1,*([0]*(1+int(f))),*row],N) for row,f in zip(values,flags)]
    c=build(N);bits=lambda v:((np.asarray(v,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    for seed in (20261028,20261029):
        sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
        for i in range(N):sim.write_across_shots(c.register_data[f'R{i}'],bits([r[i] for r in values]))
        sim.write_across_shots(c.register_data['extra'],flags[None,:].astype(bool));sim.do(c)
        for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'R{i}']),bits([r[i] for r in expected])),(N,i)
        assert np.array_equal(sim.read_across_shots(c.register_data['extra']),flags[None,:].astype(bool))
        assert not sim.read_phase_flipped_across_shots().any(),(N,'phase')
        scratch=km.array([km.q(i) for i in range(12*N+1,c.num_qubits)])
        assert not sim.read_across_shots(scratch).any(),(N,'scratch')
        sim.do(c)
        for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'R{i}']),bits([r[i] for r in values])),(N,i,'twice')
        assert not sim.read_phase_flipped_across_shots().any(),(N,'twice_phase')
        assert not sim.read_across_shots(scratch).any(),(N,'twice_scratch')
    report=dict(status='PASS_COMPONENT',N=N,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,
                scope='coherent quantum selector, both shifts, independent rational-series model, phase, scratch, and involution')
    print(json.dumps(report),flush=True);return report


if __name__=='__main__':
    rows=[]
    for N in (6,14,34):
        rows.append(check(N));Path('results/variable_shift_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
