"""Classical check of a monic Euclid state occupying exactly N field words.

At stage k>=1, retain two monic cofactors Bprev (degree k-1), B (degree k),
the band R[k:N-k], and inverse(p), where p is the previous remainder's lead.
The leading monic coefficients are implicit. The invariant is
    P*B + R*Bprev = p*x**N.
This diagnostic rejects degree gaps; it is not a submitted decoder.
"""
import json
import random
from pathlib import Path

from compact_euclid import divide, inv, mul, product, xor, triangular_products


def first(S):
    a = S[-1]
    if not a:
        raise ValueError('zero first pivot')
    pinv = inv(a)
    v = mul(S[-2], pinv)
    C = [S[j-1] ^ mul(v, S[j]) for j in range(1, len(S)-1)]
    return [1], [v, 1], C, pinv


def first_backward(Bprev, B, C, pinv):
    assert Bprev == [1] and B[-1] == 1
    a = inv(pinv)
    S = [0]*len(C) + [mul(B[0], a), a]
    for j in reversed(range(1, len(S)-1)):
        S[j-1] = C[j-1] ^ mul(B[0], S[j])
    return S


def step(k, Bprev, B, C, pinv):
    r, s = C[-1], C[-2]
    if not r:
        raise ValueError('degree gap or completion')
    newpinv = inv(r)
    beta = mul(r, pinv)
    v = mul(s, newpinv)
    alpha = v ^ B[-2]
    Bnew = xor(product([alpha, 1], B), [mul(beta, x) for x in Bprev])
    assert Bnew[-1] == 1 and Bnew[-2] == v
    out = [0]*(len(C)-2)
    for l in reversed(range(len(out))):
        degree = 2*k+1+l
        z = 0
        for i, co in enumerate(Bnew):
            j = degree-i-k
            if 0 <= j < len(C):
                z ^= mul(co, C[j])
        for h in range(1, min(k, len(out)-1-l)+1):
            z ^= mul(B[k-h], out[l+h])
        out[l] = z
    return list(B), Bnew, out, newpinv


def backward(k, B, Bnew, Cnew, newpinv):
    alpha = Bnew[-2] ^ B[-2]
    beta = Bnew[k-1] ^ (B[k-2] if k > 1 else 0) ^ mul(alpha, B[-2])
    if not beta:
        raise ValueError('singular inverse cofactor pivot')
    Bprev = [mul(x, inv(beta)) for x in xor(Bnew, product([alpha, 1], B))[:k]]
    assert Bprev[-1] == 1
    r = inv(newpinv)
    C = [0]*len(Cnew) + [mul(Bnew[-2], r), r]
    for l in reversed(range(len(Cnew))):
        degree = 2*k+1+l
        z = 0
        for i, co in enumerate(B):
            j = degree-i-(k+1)
            if 0 <= j < len(Cnew):
                z ^= mul(co, Cnew[j])
        for h in range(1, min(k+1, len(C)-1-l)+1):
            z ^= mul(Bnew[k+1-h], C[l+h])
        C[l] = z
    return Bprev, list(B), C, mul(beta, newpinv)


def monic(poly):
    iv = inv(poly[-1])
    return [mul(iv, x) for x in poly]


def estimate(N):
    rows = []
    for k in range(1, N//2):
        remain = N-2*k-2
        off = triangular_products(remain, k+1) + triangular_products(remain, k)
        boundary = sum(1 for l in range(remain) for h in [remain-l, remain+1-l] if h <= k+1)
        # Monic triangular maps have no per-word diagonal scaling. Updating
        # Bprev uses two length-k products; this conservative count does not
        # assume a fresh clean destination for their Toom combination.
        cofactor = 2*k
        scalar = 2  # beta and v, excludes inversion and phase cleanup.
        rows.append(dict(k=k, off_diagonal=off, boundary=boundary, cofactor=cofactor, scalar=scalar))
    count = N-2 + sum(sum(r[key] for key in ['off_diagonal','boundary','cofactor','scalar']) for r in rows)
    return dict(scope='serial generic-path arithmetic estimate, excludes inverse cleanup, controls, gap handling and scan', N=N,
                storage_field_words=N, field_products=count, T_field_products=42*count, rows=rows)


def test():
    rng = random.Random(919315)
    cases = []
    rejected = 0
    for N, trials in [(4,100),(6,100),(10,100),(20,30),(36,20),(70,20)]:
        good = 0
        while good < trials:
            S = [rng.randrange(4096) for _ in range(N)]
            P = [0]*N+[1]
            R = S[:]
            fp, f = [], [1]
            states = []
            try:
                state = first(S)
                for k in range(1,N//2+1):
                    q, nextR = divide(P,R)
                    if len(q) != 2:
                        raise ValueError('degree gap')
                    fp, f = f, xor(fp, product(q,f))
                    P, R = R, nextR
                    Bprev, B, C, pinv = state
                    assert Bprev == monic(fp) and B == monic(f), (N,k,'cofactor')
                    scale = inv(f[-1])
                    padded = R + [0]*N
                    assert C == [mul(scale,x) for x in padded[k:N-k]], (N,k,'remainder')
                    assert pinv == inv(mul(P[-1],inv(fp[-1]))), (N,k,'norm')
                    assert (len(Bprev)-1)+(len(B)-1)+len(C)+1 == N
                    states.append(state)
                    if k < N//2:
                        state = step(k,*state)
                        assert backward(k,*state) == states[-1], (N,k,'inverse')
                for k in reversed(range(1,N//2)):
                    state = backward(k,*state)
                    assert state == states[k-1]
                assert first_backward(*state) == S
            except ValueError:
                rejected += 1
                continue
            good += 1
        cases.append(dict(N=N,trials=trials))
    out = dict(status='PASS',scope='monic compressed Euclid, full remainder/cofactor comparisons and exact reverse',
               cases=cases,rejected_degree_gaps=rejected)
    Path('results/monic_compact_euclid_check.json').write_text(json.dumps(out,indent=2)+'\n')
    profile = estimate(70)
    Path('results/monic_compact_euclid_resources.json').write_text(json.dumps(profile,indent=2)+'\n')
    print(json.dumps(out),flush=True)
    print(json.dumps({k:v for k,v in profile.items() if k!='rows'}),flush=True)


if __name__ == '__main__':
    test()
