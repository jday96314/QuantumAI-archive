"""Reversible modular MAC using fast products and triangular division.

Reuses the attributed Apache-2.0 Gu term emitter and round-4 unit-series
middle products. These are polynomial arithmetic helpers, not the locator
recurrence. This avoids rejecting F solely on a schoolbook implementation.
"""
import json
import subprocess
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from mobius_series import BoundedArithmetic
from published_dqi import isomorphism
from oracle import product,divide

ROOT=Path(__file__).resolve().parent


def modular_mac(m,E,math,a,c,sigma,out):
    d=len(sigma)
    with m.opened(E,'poly_product',12*(2*d-1)) as flat:
        w=[flat[12*i:12*i+12] for i in range(2*d-1)]
        conv=m.polyprod_rt(a,c,[[word] for word in w])
        high=w[d:];rev=high[::-1]
        coeff=[()]+list(reversed(sigma))
        def division(inverse=False):
            if inverse:
                if high:E.emit_mac(m.polyprod_rt(sigma,high,[[w[i]] if i<d else [] for i in range(2*d-2)]))
                math.triangular(coeff,1,rev)
            else:
                math.triangular(coeff,1,rev,True)
                if high:E.emit_mac(m.polyprod_rt(sigma,high,[[w[i]] if i<d else [] for i in range(2*d-2)]))
        E.emit_mac(conv);division()
        for x,y in zip(w,out):
            for xx,yy in zip(x,y):E.cx(xx,yy)
        division(True);E.emit_mac(conv)


def f_lifetime(m,E,math,omega,psi,sigma):
    d=len(sigma)
    with m.opened(E,'old_omega',12*d) as flat:
        w=[flat[12*i:12*i+12] for i in range(d)]
        def derivative_product():
            with m.opened(E,'formal_derivative',12*d) as df:
                deriv=[df[12*i:12*i+12] for i in range(d)]
                def copy():
                    for i in range(1,d,2):
                        for x,y in zip(sigma[i],deriv[i-1]):E.cx(x,y)
                    if d%2:E.x(deriv[-1][0])
                copy();modular_mac(m,E,math,omega,deriv,sigma,w);copy()
        modular_mac(m,E,math,omega,psi,sigma,w)
        E.swaps(None,m.word_pairs(omega,w));derivative_product()
        # A phase scan could be inserted here. Retain the old-omega word
        # bank regardless of whether the reciprocal identity holds.
        derivative_product();E.swaps(None,m.word_pairs(omega,w))
        modular_mac(m,E,math,omega,psi,sigma,w)


def profile(d):
    m=load();b=m.Builder();arrays=[[b.qubits(12) for _ in range(d)] for _ in range(4)]
    class Count(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=Count(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.known_zero.add(())
    modular_mac(m,E,math,*arrays)
    macT=E.T+b._stream_totals()[0];macQ=b.peak_qubits
    bf=m.Builder();omega,psi,sig=[[bf.qubits(12) for _ in range(d)] for _ in range(3)]
    Ef=Count(bf,*bf.qubits(3),m.Alloc(bf));mf=BoundedArithmetic(m,Ef);mf.known_zero.add(())
    f_lifetime(m,Ef,mf,omega,psi,sig)
    return dict(degree=d,T=macT,Q=macQ,
                public_Q=48*d,scratch_Q=macQ-48*d,
                full_F_overwrite_lifetime_T=Ef.T+bf._stream_totals()[0],F_lifetime_Q=bf.peak_qubits,
                out_of_place_F_lifetime_T=2*macT,
                scope='exact term-emitter count; includes forward/reverse convolution, triangular division and remainder cross product')


def build_and_test(d,shots=512,kind='multiply'):
    m=load();b=m.Builder();a,c,sigma=[[b.qubits(12) for _ in range(d)] for _ in range(3)]
    out=[b.qubits(12) for _ in range(d)] if kind=='multiply' else None
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.known_zero.add(())
    if kind=='multiply':modular_mac(m,E,math,a,c,sigma,out)
    else:f_lifetime(m,E,math,a,c,sigma)
    name=f'fast_modular_d{d}' if kind=='multiply' else f'fast_F_lifetime_d{d}'
    prefixes=['a','c','sigma']+(['out'] if out is not None else [])
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tape=ROOT/'parameters'/f'{name}.npy';np.save(tape,rows)
    pb=km.CircuitBuilder()
    for prefix in prefixes:
        for i in range(d):pb.create_quantum_register(12,f'{prefix}{i}')
    pre=ROOT/'parameters'/f'{name}.pre.kmx';pb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters/empty';empty.write_text('');path=ROOT/'candidates'/f'{name}.kmx'
    subprocess.run([str(ROOT/'../../4/codex/translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(path),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    circ=km.Circuit.from_file(str(path),format='kmx');circ.to_file(str(path.with_suffix('.kmb')),format='kmb');path.unlink()
    rng=np.random.default_rng(541739+d);inputs={k:rng.integers(0,4096,(d,shots)) for k in prefixes}
    expected={k:v.copy() for k,v in inputs.items()}
    if kind=='multiply':
        for s in range(shots):
            rem=divide(product(inputs['a'][:,s],inputs['c'][:,s]),[*inputs['sigma'][:,s],1])[1]
            for i,x in enumerate(rem):expected['out'][i,s]^=x
    fw,bk,transform=isomorphism(m);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=shots,seed=99103);sim.use_same_registers_as(circ)
    for prefix,arr in inputs.items():
        for i,v in enumerate(arr):sim.write_across_shots(circ.register_data[f'{prefix}{i}'],bits(np.array([transform(int(x),fw) for x in v])))
    sim.do(circ)
    for prefix,arr in expected.items():
        for i,v in enumerate(arr):assert np.array_equal(sim.read_across_shots(circ.register_data[f'{prefix}{i}']),bits(np.array([transform(int(x),fw) for x in v]))),(d,prefix,i)
    public=set(q for reg in circ.register_data.values() for q in reg)
    scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any(),(d,'scratch')
    assert not sim.read_phase_flipped_across_shots().any(),(d,'phase')
    pr=profile(d)
    assert pr['T' if kind=='multiply' else 'full_F_overwrite_lifetime_T']==circ.operation_statistics().num_toffoli_like
    assert pr['Q' if kind=='multiply' else 'F_lifetime_Q']==circ.num_qubits
    return dict(status='PASS',kind=kind,degree=d,shots=shots,Q=circ.num_qubits,T=circ.operation_statistics().num_toffoli_like,
                scope='arbitrary monic quantum moduli and arbitrary incoming target; exact product/remainder, source restoration, clean scratch and phase')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--reuse-native-mac',action='store_true');args=p.parse_args()
    if args.reuse_native_mac:tests=json.loads(Path('results/fast_modular.json').read_text())['native_checks']
    else:tests=[build_and_test(d,256 if d==35 else 512) for d in [2,3,35]]
    tests=[r for r in tests if r.get('kind','multiply')=='multiply']+[build_and_test(3,512,'F_lifetime')]
    profiles=[profile(d) for d in [2,3,17,35]]
    out=dict(status='PASS',native_checks=tests,profiles=profiles,
             limitation='F needs 2 MACs per in-place direction, retained old-omega work and derivative construction; compact prior -> omega conversion is still extra')
    Path('results/fast_modular.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
