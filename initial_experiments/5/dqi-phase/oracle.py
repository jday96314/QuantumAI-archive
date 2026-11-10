"""Independent polynomial oracle in the HOST polynomial basis (0x1053).

No locator recurrence is implemented. Random locators are built from known
errors; general squarefree polynomials use ordinary classical polynomial EEA.
The EEA here is a correctness oracle, never a reversible gate-cost claim.
"""
import argparse
import json
import random
from pathlib import Path
import numpy as np

MODULUS = 0x1053
ORDER = 4095


def slowmul(a, b):
    out = 0
    while b:
        if b & 1: out ^= a
        b >>= 1
        a <<= 1
        if a & 4096: a ^= MODULUS
    return out


EXP = [1]
for _ in range(1, ORDER): EXP.append(slowmul(EXP[-1], 2))
assert len(set(EXP)) == ORDER
LOG = [0] * 4096
for _i, _v in enumerate(EXP): LOG[_v] = _i
EXP = np.array(EXP, dtype=np.int64)
LOG = np.array(LOG, dtype=np.int64)


def mul(a, b):
    return int(EXP[(int(LOG[a]) + int(LOG[b])) % ORDER]) if a and b else 0


def vmul(a, b):
    a, b = np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)
    return np.where((a != 0) & (b != 0), EXP[(LOG[a] + LOG[b]) % ORDER], 0)


def power(a, k):
    if k == 0: return 1
    return int(EXP[int(LOG[a]) * k % ORDER]) if a else 0


def inv(a): return power(a, 4094)


def trim(p):
    p = list(p)
    while p and p[-1] == 0: p.pop()
    return p


def add(a, b):
    c = [0] * max(len(a), len(b))
    for i, x in enumerate(a): c[i] ^= int(x)
    for i, x in enumerate(b): c[i] ^= int(x)
    return trim(c)


def product(a, b):
    if len(a) == 0 or len(b) == 0: return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] ^= mul(int(x), int(y))
    return trim(c)


def divide(a, b):
    a, b = trim(a), trim(b)
    if not b: raise ZeroDivisionError('zero polynomial')
    q = [0] * max(0, len(a) - len(b) + 1)
    vi = inv(b[-1])
    while len(a) >= len(b):
        k, v = len(a) - len(b), mul(a[-1], vi)
        q[k] = v
        for j, x in enumerate(b): a[j+k] ^= mul(v, x)
        a = trim(a)
    return trim(q), a


def inverse_mod(a, modulus):
    if len(trim(modulus)) <= 1: return []
    r, s, u, v = list(modulus), trim(a), [], [1]
    while s:
        q, rem = divide(r, s)
        r, s, u, v = s, rem, v, add(u, product(q, v))
    if len(r) != 1: raise ValueError('nonunit / repeated-root locator')
    return divide([mul(c, inv(r[0])) for c in u], modulus)[1]


def derivative(p): return trim([p[i] if i & 1 else 0 for i in range(1, len(p))])


def evaluate(p, x):
    v = np.zeros_like(np.asarray(x, dtype=np.int64))
    for c in reversed(p): v = vmul(v, x) ^ int(c)
    return int(v) if v.ndim == 0 else v


def trace(x):
    out = 0
    for _ in range(12): out ^= x; x = mul(x, x)
    assert out in (0, 1)
    return out


TRACE_MASK = sum(trace(1 << i) << i for i in range(12))
TRACE_VALUE = np.array([(x & TRACE_MASK).bit_count() & 1 for x in range(4096)])
DUAL_TO_MASK = [sum(trace(mul(u, 1 << i)) << i for i in range(12)) for u in range(4096)]
assert len(set(DUAL_TO_MASK)) == 4096
MASK_TO_DUAL = np.argsort(DUAL_TO_MASK)


def matrix_inverse(a):
    a = [list(map(int, row)) + [int(i == j) for j in range(len(a))] for i, row in enumerate(a)]
    n = len(a)
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        s = inv(a[j][j]); a[j] = [mul(s, v) for v in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                c = a[i][j]
                a[i] = [v ^ mul(c, w) for v, w in zip(a[i], a[j])]
    return np.array([r[n:] for r in a], dtype=np.int64)


def interpolation_inverse(n=70):
    return matrix_inverse([[power(int(EXP[j]), i+1) for i in range(n)] for j in range(n)])


def mask_residual(masks, inverse):
    n = len(inverse)
    mu = MASK_TO_DUAL[np.asarray(masks)]
    u = np.bitwise_xor.reduce(vmul(inverse, mu[:n][None, :]), axis=1)
    residual = mu ^ evaluate([0, *u], EXP[:len(masks)])
    assert not residual[:n].any()
    return u, residual


def known_error_polynomials(positions, values):
    sigma = [1]
    for j in positions: sigma = product(sigma, [1, int(EXP[j])])
    syndromes = [0] * 70
    for j, e in zip(positions, values):
        for i in range(70): syndromes[i] ^= mul(e, power(int(EXP[j]), i+1))
    omega = product(sigma, syndromes)[:70]
    omega = trim(omega)
    assert len(omega) <= len(positions)
    return sigma, omega, syndromes


def run(trials):
    import kickmix as km
    from import_published import load
    from published_dqi import isomorphism
    rng = random.Random(523091)
    m = load(); forward, back, transform = isomorphism(m)
    host = km.GF2Field(12)
    assert host.modulus == MODULUS and host.primitive_element == 2
    for a in range(4096):
        assert host.mul(a, 2) == slowmul(a, 2)
        assert transform(transform(a, forward), back) == a
    for _ in range(4096):
        a, b = rng.randrange(4096), rng.randrange(4096)
        assert mul(a, b) == host.mul(a, b) == slowmul(a, b)
        assert transform(m.gmul(transform(a, forward), transform(b, forward)), back) == mul(a, b)
    vi = interpolation_inverse()
    checks = roots = masks_checked = 0; skip_hist = {}; compact_checks = 0
    for degree in range(36):
        for trial in range(trials):
            pos = rng.sample(range(4095), degree)
            vals = [rng.randrange(1, 4096) for _ in pos]
            sigma, omega, syndromes = known_error_polynomials(pos, vals)
            psi = inverse_mod(derivative(sigma), sigma)
            F = divide(product(omega, psi), sigma)[1]
            if degree:
                assert divide(product(psi, derivative(sigma)), sigma)[1] == [1]
            beta = EXP[(-np.arange(4095)) % 4095]
            rv = evaluate(sigma, beta)
            assert set(np.flatnonzero(rv == 0)) == set(pos)
            fv, ov, dv, pv = [evaluate(p, beta) for p in (F, omega, derivative(sigma), psi)]
            for j, e in zip(pos, vals):
                assert fv[j] == e == mul(int(ov[j]), inv(int(dv[j])))
                assert mul(int(pv[j]), int(dv[j])) == 1
                roots += 1
            masks = np.array([rng.randrange(4096) for _ in range(4095)])
            if trial == 0:
                masks[:] = 0 if degree == 0 else masks
                if degree == 1: masks[:] = 4095
            want = sum((e & int(masks[j])).bit_count() for j, e in zip(pos, vals)) & 1
            got = sum((int(fv[j]) & int(masks[j])).bit_count() for j in np.flatnonzero(rv == 0)) & 1
            assert want == got
            u, residual = mask_residual(masks, vi)
            zlin = sum(trace(mul(int(a), s)) for a, s in zip(u, syndromes)) & 1
            residual_phase = sum(trace(mul(int(residual[j]), int(fv[j]))) for j in np.flatnonzero(rv == 0)) & 1
            assert want == zlin ^ residual_phase
            skips = int(np.count_nonzero(residual == 0)); skip_hist[skips] = skip_hist.get(skips, 0) + 1
            # The compact cofactor is algebraically recovered from omega;
            # no recurrence is used to establish this scan identity.
            if degree:
                xn = divide([0] * 70 + [1], sigma)[1]
                scaled_bp = divide(product(xn, inverse_mod(omega, sigma)), sigma)[1]
                for j, e in zip(pos, vals):
                    a = evaluate(scaled_bp, int(beta[j])); d = int(dv[j])
                    assert mul(power(int(beta[j]), 70), inv(mul(a, d))) == e
                compact_checks += 1
            checks += 1; masks_checked += 4095
    # Non-split squarefree polynomials: reciprocal identity remains valid.
    general = nonsquarefree = 0
    for degree in range(1, 36):
        p = [rng.randrange(4096) for _ in range(degree)] + [1]
        try: psi = inverse_mod(derivative(p), p)
        except ValueError: nonsquarefree += 1; continue
        assert divide(product(psi, derivative(p)), p)[1] == [1]
        general += 1
    repeated = product([1, 2], [1, 2])
    try: inverse_mod(derivative(repeated), repeated)
    except ValueError: pass
    else: raise AssertionError('repeated-root inverse accepted')
    zero_root = product([0, 1], [1, 2])
    assert divide(product(inverse_mod(derivative(zero_root), zero_root), derivative(zero_root)), zero_root)[1] == [1]
    out = dict(status='PASS', scope='independent host-basis polynomial and mask oracle; not native decoder acceptance', seed=523091,
        trials_per_degree=trials, split_cases=checks, root_checks=roots, arbitrary_mask_words=masks_checked,
        general_squarefree_cases=general, general_nonunit_cases=nonsquarefree, compact_inverse_product_cases=compact_checks,
        mask_skip_histogram=skip_hist, exceptional_cases=['degree zero', 'zero-root squarefree locator', 'repeated-root rejection', 'all-zero mask', 'all-one mask'],
        host_modulus=hex(MODULUS), position='j -> alpha^(-j), alpha=2',
        identity='e_j=omega(beta)/sigma_prime(beta)=F(beta)=beta^70/(scaled_Bprev(beta)*Bprime(beta))',
        trace_dual_columns=[int(MASK_TO_DUAL[1<<i]) for i in range(12)], upstream_modulus=hex(m.POLY))
    Path('results/oracle.json').write_text(json.dumps(out, indent=2)+'\n'); print(json.dumps(out), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--trials', type=int, default=4); a = p.parse_args(); run(a.trials)
