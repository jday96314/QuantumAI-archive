"""Native variable-degree generalized inverse and arbitrary paired inverse."""
import json,time,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
import padded_polynomial_inverse as pad
import check_blocked_cfraction as ref

HERE=Path(__file__).resolve().parent

def build(n,both=False):
    m=pad.base.load();m.inv_word=pad.base.bind(m);b=m.Builder()
    F=[b.qubits(12) for _ in range(n+1)];G=[b.qubits(12) for _ in range(n)];flags=b.qubits(n)
    E=m.Em(b,*b.qubits(3),m.Alloc(b));pad.apply(m,E,F,G,flags)
    if both:pad.apply(m,E,F,G,flags,True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tag=f'padded_poly_inverse_n{n}'+('_both' if both else '')
    tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for name in ('F','G'):
        for i in range(n+int(name=='F')):kb.create_quantum_register(12,f'{name}{i}')
    kb.create_quantum_register(n,'flags');pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(n,shots=256):
    start=time.time();m=pad.base.load();rng=np.random.default_rng(2026110300+n);rows=[]
    for both in (False,True):
        fs=[];gs=[];flags=[];want=[]
        for j in range(shots):
            d=j%(n+1);z=n-d
            if both:
                f=rng.integers(0,4096,n+1).tolist();g=rng.integers(0,4096,n).tolist()
                ff=f;gg=g;fl=rng.integers(0,2,n).astype(bool).tolist()
            else:
                while True:
                    f=rng.integers(0,4096,d).tolist()+[1];g=rng.integers(0,4096,d).tolist()
                    if not d or pad.base.model(m,f,g)[2]:break
                f=[0]*z+f;g=[0]*z+g;ff,gg,_=pad.base.model(m,f,g,degree=d)
                fl=[i>=d for i in range(n)]
            fs.append(f);gs.append(g);flags.append(fl);want.append((ff,gg))
        c=build(n,both);bits=lambda a:((np.asarray(a,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
        for seed in (20261103,20261104):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for name,values in (('F',fs),('G',gs)):
                for i in range(n+int(name=='F')):sim.write_across_shots(c.register_data[f'{name}{i}'],bits([r[i] for r in values]))
            sim.write_across_shots(c.register_data['flags'],np.asarray(flags,dtype=bool).T);sim.do(c)
            for k,name in enumerate(('F','G')):
                for i in range(n+int(name=='F')):assert np.array_equal(sim.read_across_shots(c.register_data[f'{name}{i}']),bits([r[k][i] for r in want])),(n,both,name,i)
            assert np.array_equal(sim.read_across_shots(c.register_data['flags']),np.asarray(flags,dtype=bool).T)
            assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(25*n+12,c.num_qubits)])).any(),(n,both,'scratch')
        rows.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots))
    result=dict(status='PASS_COMPONENT',n=n,circuits=rows,seconds=time.time()-start,
                scope='forward all degrees with common x padding and regular active pivots; paired inverse arbitrary polynomials/flags; phase and scratch')
    print(json.dumps(result),flush=True);return result

if __name__=='__main__':
    rows=[]
    for n in (3,7,17):
        rows.append(test(n));(HERE/'results/padded_poly_inverse_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
