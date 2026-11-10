"""Independent integer model for fused rounded-square inverse approximation."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
TABLE = np.load(PROJECT / 'initial_experiments/4/codex/candidates/margin_curve8_t5.9e-08.npy')
MAX_X = 3133452


def rounded_product(a, coordinate, width, exponent):
    result = np.zeros_like(coordinate)
    bit = lambda j: 0 if j < 0 else (coordinate >> min(j, width - 1)) & 1
    for g in range((width + 1) // 2):
        digit = bit(2*g - 1) + bit(2*g) - 2*bit(2*g + 1)
        for mag in (1, 2):
            e = exponent + 2*g + (mag == 2)
            v = a << e if e >= 0 else (a >> -e) + ((a >> (-e-1)) & 1)
            result += np.where(abs(digit) == mag, np.sign(digit) * v, 0)
    return result


def evaluate(xs, fractional=2, q=5, square_round=16, detail=False):
    xs = np.asarray(xs, dtype=np.int64)
    k = np.zeros_like(xs)
    for j in range(21):
        k += xs < (1 << (21-j))
    m = xs << k
    s = k // 2
    address = ((m >> 13) & 255) + ((k & 1) << 8)
    c0, c1, c2 = TABLE[address].T
    z = (m & 16383) - 8192
    u = (z + (1 << (q-1))) >> q
    square_shift = square_round - 2*q
    h = ((u*u + (1 << (square_shift-1))) >> square_shift)
    linear = rounded_product(c1, z & 16383, 14, fractional-12)
    curvature = rounded_product(c2, h, int(((256**2 + 32) >> 6)).bit_length()+1,
                                fractional+square_round-23)
    p = (c0 << fractional) - linear + curvature
    phase = p.astype(np.uint64) << (s + 27 - fractional).astype(np.uint64)
    if detail:
        return phase, dict(k=k,m=m,s=s,address=address,z=z,u=u,h=h,p=p,
                           linear=linear,curvature=curvature)
    return phase


def certify(fractional, q=5, square_round=16):
    target_path = PROJECT / 'initial_experiments/5/inv-phase/reference/exact_target.npz'
    with np.load(target_path) as d:
        target, tol = d['target'], d['tol']
    margin = 1 << 63
    worst_ratio = 0
    worst_x = 0
    for first in range(1, MAX_X+1, 131072):
        xs = np.arange(first, min(MAX_X+1, first+131072), dtype=np.int64)
        phase = evaluate(xs, fractional, q, square_round)
        error = np.abs((phase-target[first-1:first-1+len(xs)]).view(np.int64))
        t = tol[first-1:first-1+len(xs)]
        ratios = error/t
        i = int(ratios.argmax())
        if ratios[i] > worst_ratio:
            worst_ratio, worst_x = float(ratios[i]), int(xs[i])
        margin = min(margin, int((t.astype(np.int64)-error).min()))
    return dict(status='MODEL_PASS' if margin > 0 else 'MODEL_FAIL', inputs=MAX_X,
                fractional=fractional,q=q,square_round=square_round,
                worst_ratio=worst_ratio,worst_x=worst_x,min_margin_units=margin,
                min_margin_turns=margin/2**64)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fractional', type=int, nargs='+', default=[0,1,2,3,4])
    args = parser.parse_args()
    reports = []
    for f in args.fractional:
        r = certify(f)
        reports.append(r)
        print(json.dumps(r), flush=True)
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/square_models.json').write_text(json.dumps(reports, indent=2)+'\n')
