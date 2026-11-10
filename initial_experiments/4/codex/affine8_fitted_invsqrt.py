"""Affine storage for finer coordinates with aggressively narrow curvature.

Research builder only. R0=C0+q*2^17-11*2^23-p*2^25, q is the eight-bit
coarse address and p normalization parity. The original phase is restored by
11/16+p/4-q/1024. Every candidate needs native accuracy validation.
"""
import argparse
import inspect
import json
from pathlib import Path

import kickmix as km
import numpy as np

import affine_invsqrt as affine
import binary_centered_invsqrt as poly
from monomial_unlookup import erase
from optimize_workspace import compact
from short_adder import add


def padded_erase(b, address, target, table, extra_outcomes=()):
    # The added zero targets carry no table bits. Their random HMR outcomes
    # therefore contribute no phase; their wires supply the monomial banks.
    n = len(address)
    needed = sum((1 << k)-k-1 for k in [n//2, n-n//2])
    pad = b.alloc_qubits(max(0, needed-len(target)))
    erase(b, address, km.array([*target, *pad]), table, extra_outcomes)
    b.free(pad)


def make_builder(raw, widths, correction=0):
    lookup = (len(raw)//2).bit_length()-1
    if lookup > 8 and np.array_equal(raw[::2], raw[1::2]):
        raw = raw[::2]
        lookup -= 1
    index = np.arange(len(raw))
    parity = index >> lookup
    coarse = (index & ((1 << lookup)-1)) >> (lookup-8)
    tab = raw.copy()
    tab[:, 0] += (coarse << 17) - (11 << 23) - (parity << 25) - correction
    assert tab.min() >= 0 and int(tab[:, 0].max()).bit_length() <= 23

    def restore(b, x, y, regs, flags, limit, inverse=False):
        # Same cheap constant 3/16 plus borrowed normalized 1/2 as the
        # selected affine circuit; only the coarse address subtraction grows.
        source = inspect.getsource(affine.restore_affine)
        source = source.replace('x[14:21],y[18:27]', 'x[13:21],y[17:27]')
        ns = dict(affine.__dict__)
        exec(source, ns)
        ns['restore_affine'](b, x, y, regs, flags, limit, inverse)
        assert correction == 0

    source = inspect.getsource(affine.build)
    assert source.count("namespace['build'](7,(27,18,8),") == 1
    source = source.replace("namespace['build'](7,(27,18,8),", "namespace['build'](8,FIT_WIDTHS,")
    source = source.replace('lookup_frac=8', 'lookup_frac=FIT_LOOKUP')
    source = source.replace('accum_bits=[18,9,None]', 'accum_bits=FIT_SPLITS')
    # Patch the copied polynomial builder before the affine source rewrites
    # it, preserving all guarded replacement counts in the production helper.
    psource = inspect.getsource(poly.build)
    old = 'erase(b,address,target,packed,extra_outcomes=delayed_meas if delay_high_bits else ())'
    assert psource.count(old) == 1
    psource = psource.replace(old, 'padded_erase(b,address,target,packed,extra_outcomes=delayed_meas if delay_high_bits else ())')
    pns = dict(poly.__dict__)
    pns['padded_erase'] = padded_erase
    exec(psource, pns)
    # inspect.getsource needs a real source file for the affine helper's copy.
    generated = Path(f'parameters/affine8_poly_source_{widths[1]}_{widths[2]}.py')
    generated.write_text(psource)
    import importlib.util
    spec = importlib.util.spec_from_file_location('affine8_poly_source', generated)
    module = importlib.util.module_from_spec(spec)
    module.__dict__.update(pns)
    exec(compile(psource, str(generated.resolve()), 'exec'), module.__dict__)
    namespace = dict(affine.__dict__)
    namespace.update(table=lambda: tab, restore_affine=restore, p=module,
                     FIT_WIDTHS=widths, FIT_LOOKUP=lookup, FIT_SPLITS=[17, widths[1]-10, None])
    exec(source, namespace)
    return namespace['build']


def check(raw, widths):
    lf = (len(raw)//2).bit_length()-1
    old = poly.table_for
    poly.table_for = lambda *args: raw
    lo, hi, worst = 2**62, -2**62, 0.
    try:
        for start in range(1, 3133453, 262144):
            xs = np.arange(start, min(3133453, start+262144))
            k = 21-np.floor(np.log2(xs)).astype(np.int64)
            y = poly.model(xs, 8, widths, rounded=True, lookup_frac=lf, booth=True) >> (k//2)
            coarse = ((xs << k) >> 13) & 255
            residual = y + (coarse << 17) - (11 << 23) - ((k & 1) << 25)
            lo, hi = min(lo, int(residual.min())), max(hi, int(residual.max()))
            error = y.astype(float)*np.exp2(k//2)/2**37*np.sqrt(xs)-1
            worst = max(worst, float(abs(error).max()))
    finally:
        poly.table_for = old
    assert 0 <= lo <= hi < 1 << 23, (lo, hi)
    assert worst < 2**-24, worst
    return dict(minimum=lo, maximum=hi, max_relative_error=worst)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--table', required=True)
    p.add_argument('--w1', type=int, required=True)
    p.add_argument('--w2', type=int, required=True)
    p.add_argument('--tag', required=True)
    a = p.parse_args()
    raw = np.loadtxt(a.table, skiprows=1, dtype=np.int64)
    widths = (27, a.w1, a.w2)
    bounds = check(raw, widths)
    builder = make_builder(raw, widths)
    rows = []
    for limit in [78, 80, 82, 84, 86, 88, 90]:
        for keep in [0, 1, 2, 3]:
            if keep >= int(raw[:, 2].max()).bit_length():
                continue
            try:
                c = compact(builder(limit=limit, keep=keep, delay=2, product=17, deferred=3), Path(f'candidates/{a.tag}_tmp.kmb'))
            except (ValueError, AssertionError):
                continue
            T = c.operation_statistics().num_toffoli_like
            estimate = c.num_qubits*(T-11+251)
            row = dict(Q=c.num_qubits, T=T, estimate=estimate, limit=limit, keep=keep)
            if not rows or estimate < min(r['estimate'] for r in rows):
                c.to_file(f'candidates/{a.tag}_best.kmb', format='kmb')
            rows.append(row)
    rows.sort(key=lambda x:x['estimate'])
    out = dict(widths=widths, bounds=bounds, rows=rows)
    Path(f'results/{a.tag}_resources.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(dict(bounds=bounds, best=rows[:3])), flush=True)


if __name__ == '__main__':
    main()
