"""Independent polynomial and quantum checks of deferred completion padding."""
import argparse
import json
import random
import subprocess
import time
from pathlib import Path

import kickmix as km
import numpy as np

from blocked_euclid_sequence import frames as block_frames
from blocked_monic_euclid import apply as block_apply
from block_euclid_ratio import series_divide,series_product
from check_compact_forney import derivative, evaluate
from compact_euclid import divide, inv, mul, product, xor
from coprime_padding import H_of, M
from deferred_euclid_padding import apply as postprocess, direct_repad,frames, padding_polynomial
from import_published import load
from mobius_series import BoundedArithmetic
from monic_compact_euclid import first, step
from monic_euclid_local import first_step

ROOT = Path(__file__).resolve().parent


def trim(p):
    p = list(p)
    while p and not p[-1]:
        p.pop()
    return p


def syndrome(N, positions, errors):
    S = [0]*N
    for position, error in zip(positions, errors):
        value = error
        for j in range(N):
            S[j] ^= value
            value = mul(value, position)
    return S


def monomial_reference(k,count,bp,b,c,pinv):
    n=2*count
    oldb,oldc=list(b),list(c)
    t=list(reversed(c))[:n]
    f=series_divide(series_product(t,list(reversed(b)),n),[inv(pinv),0]+series_product(t,list(reversed(bp)),max(0,n-2)),n)
    flags=[]
    for j in range(count):
        zero=not f[0]
        flags.append(zero)
        if zero:f[0]=1
        beta=f[0]
        u=series_divide([beta],f,len(f))
        alpha=u[1]
        f=u[2:]
        if zero:
            if alpha:raise ValueError('non-generic degree gap')
            bp,b=[0]+bp,[0]+b
        else:
            bp,b=b,xor(product([alpha,1],b),[mul(beta,x) for x in bp])
        pinv=mul(pinv,inv(beta))
    numerator=series_product(list(reversed(b)),list(reversed(oldc)),len(oldc))
    c=list(reversed(series_divide(numerator[2*count:],list(reversed(oldb)),len(oldc)-2*count)))
    return (bp,b,c,pinv),flags


def reference(S, weight, width):
    state = first(S)
    for k in range(1, weight):
        state = step(k, *state)
    bp0, b0, c0, p0 = state
    assert not any(c0)
    state = first(S)
    flags = [False]
    k = 1
    while k < len(S)//2:
        count = min(width, len(S)//2-k)
        state, ff = monomial_reference(k, count, *state)
        flags.extend(ff)
        k += count
    if flags != [False]*weight+[True]*(len(S)//2-weight):
        raise ValueError(('non-generic completion flags', weight, flags))
    bp, b, c, pinv = state
    # Independently divide out the core's common powers of x.
    for k in reversed(range(weight, len(S)//2)):
        assert bp[0]==0 and b[0]==0
        bp,b=bp[1:],b[1:]
    assert bp == bp0 and b == b0 and pinv == p0
    P = padding_polynomial(weight, len(S)//2)
    padded_bp, padded_b = product(P, bp), product(P, b)
    H0, Hp = H_of(bp, b), H_of(padded_bp, padded_b)
    assert Hp == product(P, H0)
    recovered, remainder = divide(Hp, P)
    assert not remainder and recovered == H0
    f, g = padded_b, trim(derivative(padded_b))
    while g:
        _, remainder = divide(f, g)
        f, g = g, trim(remainder)
    assert len(f) == 1, 'padded locator is not squarefree'
    return (padded_bp, padded_b, pinv), flags, H0


def classical(N, width, trials):
    rng = random.Random(901725)
    began = time.time()
    cases = roots = gaps = 0
    for weight in range(1, N//2+1):
        for _ in range(trials):
            positions = rng.sample(range(1, 4096), weight)
            errors = [rng.randrange(1, 4096) for _ in positions]
            S = syndrome(N, positions, errors)
            try:
                (bp, b, pinv), flags, H = reference(S, weight, width)
            except (ValueError, ZeroDivisionError):
                gaps += 1
                continue
            for position, error in zip(positions, errors):
                root = inv(position)
                assert evaluate(b, root) == 0
                h = evaluate(H, root)
                got = mul(mul(inv(pinv), M.gpow(root, N-1)), M.gpow(inv(h), 2))
                assert got == error
                roots += 1
            cases += 1
    out = dict(status='PASS', scope='classical core monomial completion padding, final squarefree common-factor padding, exact H division and root magnitudes; generic paths only',
               N=N, block=width, trials_per_weight=trials, cases=cases, roots=roots,
               skipped_degree_gaps=gaps, seconds=time.time()-began)
    tag=f'deferred_padding_n{N}_b{width}'+(f'_trials{trials}' if trials!=100 else '')
    Path(f'results/{tag}_classical.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out), flush=True)


def build(N, width, which,direct=False,thin=False,suffix=False,serial_after=None):
    m = load()
    b = m.Builder()
    words = [b.qubits(12) for _ in range(N)]
    flags = b.qubits(N//2)
    E = m.Em(b, *b.qubits(3), m.Alloc(b))
    math = BoundedArithmetic(m, E)
    stages, _ = block_frames(words, width)
    if which in ['forward', 'both']:
        with m.eq_test(E, words[-1], 0) as top:
            E.cx(top, flags[0])
        E.cx(flags[0], words[-1][0])
        first_step(m, E, math, words)
        for k, count, bp, bb, c, p in stages:
            if serial_after is not None and k>=serial_after:
                from monomial_serial_euclid import apply as serial
                for j in range(count):bp,bb,c,p=serial(m,E,math,bp,bb,c,p,flags[k+j])
            else:block_apply(m, E, math, bp, bb, c, p, count, flags=flags[k:k+count],monomial=True,thin=thin,suffix=suffix)
        (direct_repad if direct else postprocess)(m, E, words, flags)
    if which in ['inverse', 'both']:
        (direct_repad if direct else postprocess)(m, E, words, flags, True)
        for k, count, bp, bb, c, p in reversed(stages):
            if serial_after is not None and k>=serial_after:
                from monomial_serial_euclid import apply as serial
                from blocked_monic_euclid import schedule
                Q=list(reversed(c))[:2*count];ff,_=schedule(bp,bb,Q)
                for j in reversed(range(count)):
                    oldprev,old,be,al=ff[j]
                    serial(m,E,math,oldprev,old,c[:len(c)-2*j],p,flags[k+j],True)
            else:block_apply(m, E, math, bp, bb, c, p, count, True, flags=flags[k:k+count],monomial=True,thin=thin,suffix=suffix)
        first_step(m, E, math, words, True)
        E.cx(flags[0], words[-1][0])
        with m.eq_test(E, words[-1], 0) as top:
            E.cx(top, flags[0])
    rows = np.asarray(b._items._rows, dtype=np.uint32).reshape(-1, 9)
    tag = f'deferred_padding'+('_direct' if direct else '')+('_thin' if thin else '')+('_suffix' if suffix else '')+(f'_serial{serial_after}' if serial_after is not None else '')+f'_n{N}_b{width}_{which}'
    tape = ROOT/f'parameters/{tag}.npy'
    np.save(tape, rows)
    kb = km.CircuitBuilder()
    for i in range(N):
        kb.create_quantum_register(12, f'S{i}')
    kb.create_quantum_register(N//2, 'zero_pivots')
    pre = ROOT/f'parameters/{tag}.pre.kmx'
    kb.finish_circuit().to_file(str(pre), format='kmx')
    empty = ROOT/'parameters/deferred_padding.empty'
    empty.write_text('')
    out = ROOT/f'candidates/{tag}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'), str(tape), str(len(rows)), '0', str(pre), str(empty), str(empty), str(out), str(b.next_out), str(empty), str(len(rows))], check=True)
    c = km.Circuit.from_file(str(out), format='kmx')
    c.to_file(str(out.with_suffix('.kmb')), format='kmb')
    return c


def quantum(N, width, shots,direct=False,thin=False,suffix=False,serial_after=None):
    rng = random.Random(718930)
    initial, expected, flagvalues = [], [], []
    _, (bpidx, bidx) = frames(list(range(N)))
    gaps = 0
    while len(initial) < shots:
        weight = 1+len(initial) % (N//2)
        positions = rng.sample(range(1, 4096), weight)
        errors = [rng.randrange(1, 4096) for _ in positions]
        S = syndrome(N, positions, errors)
        try:
            (bp, b, pinv), flags, _ = reference(S, weight, width)
        except (ValueError, ZeroDivisionError):
            gaps += 1
            continue
        target = S.copy()
        for j, value in zip(bpidx, bp[:-1]):
            target[j] = value
        for j, value in zip(bidx, b[:-1]):
            target[j] = value
        target[-1] = pinv
        initial.append(S)
        expected.append(target)
        flagvalues.append(flags)
    initial = np.array(initial).T
    expected = np.array(expected).T
    flagvalues = np.array(flagvalues, dtype=bool).T
    bits = lambda v: ((v[None, :] >> np.arange(12)[:, None]) & 1).astype(bool)
    records = []
    for which in ['forward', 'inverse', 'both']:
        c = build(N, width, which,direct,thin,suffix,serial_after)
        sim = km.Simulator(batch_size=shots, seed=71937)
        sim.use_same_registers_as(c)
        start = expected if which == 'inverse' else initial
        want = expected if which == 'forward' else initial
        if which == 'both':
            # Fresh arbitrary inputs, including zeros and deliberately missing
            # top pivots, check exact reversal independently of generic decoding.
            start = np.array([[rng.randrange(4096) for _ in range(shots)] for _ in range(N)])
            start[:, 0] = 0
            start[-1, 1::8] = 0
            want = start
        for j, values in enumerate(start):
            sim.write_across_shots(c.register_data[f'S{j}'], bits(values))
        if which == 'inverse':
            sim.write_across_shots(c.register_data['zero_pivots'], flagvalues)
        sim.do(c)
        for j, values in enumerate(want):
            assert np.array_equal(sim.read_across_shots(c.register_data[f'S{j}']), bits(values)), (which, j)
        wantedflags = flagvalues if which == 'forward' else np.zeros_like(flagvalues)
        assert np.array_equal(sim.read_across_shots(c.register_data['zero_pivots']), wantedflags)
        public = {q for reg in c.register_data.values() for q in reg}
        scratch = km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(), (which, 'scratch')
        assert not sim.read_phase_flipped_across_shots().any(), (which, 'phase')
        records.append(dict(which=which, Q=c.num_qubits, T=c.operation_statistics().num_toffoli_like))
    record = dict(status='PASS', scope='full armed transform plus deferred common-factor padding; independent polynomial values, inverse, arbitrary-input composition, flags, scratch and phase; no Wronskian or decoder phase',
                  N=N, block=width, direct=direct,thin=thin,suffix=suffix,serial_after=serial_after,shots=shots, skipped_degree_gaps=gaps, circuits=records)
    tag='deferred_padding'+('_direct' if direct else '')+('_thin' if thin else '')+('_suffix' if suffix else '')+(f'_serial{serial_after}' if serial_after is not None else '')+f'_n{N}_b{width}'
    Path(f'results/{tag}_quantum.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--n', type=int, default=22)
    p.add_argument('--block', type=int, default=4)
    p.add_argument('--trials', type=int, default=100)
    p.add_argument('--shots', type=int, default=512)
    p.add_argument('--classical-only', action='store_true')
    p.add_argument('--direct',action='store_true')
    p.add_argument('--thin',action='store_true')
    p.add_argument('--suffix',action='store_true')
    p.add_argument('--serial-after',type=int)
    a = p.parse_args()
    classical(a.n, a.block, a.trials)
    if not a.classical_only:
        quantum(a.n, a.block, a.shots,a.direct,a.thin,a.suffix,a.serial_after)
