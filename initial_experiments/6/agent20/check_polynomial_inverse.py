"""Native continuation-reversal inverse with arbitrary and zero-pivot states."""
import json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
import polynomial_inverse_reversal as inv
import check_blocked_cfraction as ref

HERE=Path(__file__).resolve().parent

def build(n):
    m=inv.load();m.inv_word=inv.bind(m);b=m.Builder()
    F=[b.qubits(12) for _ in range(n)];G=[b.qubits(12) for _ in range(n)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));inv.apply(m,E,F,G)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tag=f'poly_inverse_reversal_d{n}';tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for name in ('F','G'):
        for i in range(n):kb.create_quantum_register(12,f'{name}{i}')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(n,shots=512):
    start=time.time();m=inv.load();rng=np.random.default_rng(2026103100+n)
    fs=rng.integers(0,4096,(shots,n)).tolist();gs=rng.integers(0,4096,(shots,n)).tolist()
    fs[0]=[0]*n;gs[0]=[0]*n;fs[1]=[1]*n;gs[1]=[1]*n
    for i in range(2,shots,4):gs[i][-1]=0
    want=[inv.model(m,f,g) for f,g in zip(fs,gs)]
    for f,g,(ff,gg,ok) in zip(fs,gs,want):
        assert inv.model(m,ff,gg)[:2]==(f,g)
        if ok:
            assert f==ff
            prod=[0]*(2*n-1)
            for i,x in enumerate(g):
                for j,y in enumerate(gg):prod[i+j]^=m.gmul(x,y)
            for i in reversed(range(n,len(prod))):
                for j,x in enumerate(f):prod[i-n+j]^=m.gmul(prod[i],x)
            assert prod[:n]==[1]+[0]*(n-1)
    c=build(n);bits=lambda a:((np.asarray(a,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    for seed in (20261031,20261101):
        sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
        for name,values in (('F',fs),('G',gs)):
            for i in range(n):sim.write_across_shots(c.register_data[f'{name}{i}'],bits([r[i] for r in values]))
        for expect in ([(f,g) for f,g,_ in want],list(zip(fs,gs))):
            sim.do(c)
            for k,name in enumerate(('F','G')):
                for i in range(n):assert np.array_equal(sim.read_across_shots(c.register_data[f'{name}{i}']),bits([r[k][i] for r in expect])),(n,name,i)
            assert not sim.read_phase_flipped_across_shots().any(),(n,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(24*n,c.num_qubits)])).any(),(n,'scratch')
    report=dict(status='PASS_COMPONENT',degree=n,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,
                generic_charts=sum(bool(ok) for _,_,ok in want),seconds=time.time()-start,
                scope='arbitrary totalized-model match, phase, scratch and involution; product-modulus identity on nonzero pivots')
    print(json.dumps(report),flush=True);return report

if __name__=='__main__':
    rows=[]
    for n in (1,2,3,7,17,35):
        rows.append(test(n));(HERE/'results/poly_inverse_reversal_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
