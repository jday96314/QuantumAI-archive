"""Native projective inverse: all padding degrees and arbitrary paired inverse."""
import argparse,json,time,hashlib,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
import projective_polynomial_inverse as projective
import check_blocked_cfraction as ref
HERE=Path(__file__).resolve().parent


def build(n,both=False):
    m=projective.base.load();m.inv_word=projective.bind(m);b=m.Builder()
    F=[b.qubits(12) for _ in range(n+1)];G=[b.qubits(12) for _ in range(n)];flags=b.qubits(n)
    E=m.Em(b,*b.qubits(3),m.Alloc(b));projective.padded_apply(m,E,F,G,flags)
    if both:projective.padded_apply(m,E,F,G,flags,True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tag=f'projective2_poly_inverse_n{n}'+('_both' if both else '')
    tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for name in ('F','G'):
        for i in range(n+int(name=='F')):kb.create_quantum_register(12,f'{name}{i}')
    kb.create_quantum_register(n,'flags');pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def test(n,shots=256):
    started=time.time();m=projective.base.load();rng=np.random.default_rng(2026102504+n);rows=[]
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
                    if not d:scalar=1;ginv=[];break
                    try:scalar=projective.gcd_constant(m,f,g)
                    except ValueError:continue
                    _,ginv,ok=projective.base.model(m,f,g)
                    if ok and scalar:break
                si2=m.gmul(m.ginv(scalar),m.ginv(scalar))
                ff=[0]*z+[m.gmul(si2,v) for v in f];gg=[0]*z+ginv
                f=[0]*z+f;g=[0]*z+g;fl=[i>=d for i in range(n)]
            fs.append(f);gs.append(g);flags.append(fl);want.append((ff,gg))
        c=build(n,both);bits=lambda a:((np.asarray(a,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
        for seed in (2026102505,2026102506):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for name,values in (('F',fs),('G',gs)):
                for i in range(n+int(name=='F')):sim.write_across_shots(c.register_data[f'{name}{i}'],bits([r[i] for r in values]))
            sim.write_across_shots(c.register_data['flags'],np.asarray(flags,dtype=bool).T);sim.do(c)
            for k,name in enumerate(('F','G')):
                for i in range(n+int(name=='F')):
                    actual=sim.read_across_shots(c.register_data[f'{name}{i}']);expected=bits([r[k][i] for r in want])
                    assert np.array_equal(actual,expected),(n,both,name,i,np.flatnonzero(np.any(actual^expected,axis=0))[:8])
            assert np.array_equal(sim.read_across_shots(c.register_data['flags']),np.asarray(flags,dtype=bool).T)
            assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(25*n+12,c.num_qubits)])).any(),(n,both,'scratch')
        path=HERE/'candidates'/f'projective2_poly_inverse_n{n}{"_both" if both else ""}.kmb'
        rows.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    result=dict(status='PASS_COMPONENT',n=n,circuits=rows,seconds=time.time()-started,
                scope='projective locator scale and true modular inverse on all degrees/regular active pivots; arbitrary paired inverse and non-unary flags; phase and scratch')
    print(json.dumps(result),flush=True);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degrees',nargs='+',type=int,default=[3,7,17]);a=p.parse_args();rows=[]
    for n in a.degrees:
        rows.append(test(n));(HERE/'results/projective2_poly_inverse_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
