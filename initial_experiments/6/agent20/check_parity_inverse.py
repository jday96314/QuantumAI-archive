"""Independent exact/native reciprocal tests, including implicit constant."""
import json,subprocess
from pathlib import Path
import kickmix as km
import numpy as np
from import_published import load
from parity_series_inverse import ParityArithmetic
import check_blocked_cfraction as ref

HERE=Path(__file__).resolve().parent

def build(N,zero1=False):
    m=load();b=m.Builder();F=[b.qubits(12) for _ in range(N)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=ParityArithmetic(m,E)
    if zero1:math.known_zero.add(tuple(F[1]))
    math.inverse(F)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tag=f'parity_inverse_n{N}'+('_zero1' if zero1 else '')
    tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'F{i}')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(N,zero1=False,shots=512):
    m=load();rng=np.random.default_rng(2026103000+N);src=rng.integers(0,4096,(shots,N)).tolist()
    src[:2]=[[0]*N,[1]*N]
    if zero1:
        for row in src:row[1]=0
    want=[[row[0],*ref.divide(m,[1],[1,*row[1:]],N)[1:]] for row in src]
    c=build(N,zero1);bits=lambda a:((np.asarray(a,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    for seed in (20261030,20261031):
        sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
        for i in range(N):sim.write_across_shots(c.register_data[f'F{i}'],bits([r[i] for r in src]))
        for expected in (want,src):
            sim.do(c)
            for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'F{i}']),bits([r[i] for r in expected])),(N,zero1,i)
            assert not sim.read_phase_flipped_across_shots().any(),(N,zero1,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(12*N,c.num_qubits)])).any(),(N,zero1,'scratch')
    row=dict(status='PASS_COMPONENT',N=N,zero1=zero1,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,
             scope='exact independent reciprocal model, arbitrary ignored constant word, phase, scratch, involution')
    print(json.dumps(row),flush=True);return row

if __name__=='__main__':
    rows=[]
    for N in (2,3,6,8,14,34):
        for zero1 in (False,True):
            rows.append(test(N,zero1));(HERE/'results/parity_inverse_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
