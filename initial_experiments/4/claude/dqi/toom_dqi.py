"""Generic decoder whose step update multiplies by q = q0 + q1 x via evaluation/interpolation.

c = q * y (deg y = n - 1) needs n + 1 products: p_inf = q1 y_{n-1} and
p_j = q(a_j) y(a_j) at n distinct points a_j. Evaluation, interpolation and
constant multiplications over GF(4096) are CNOT networks; each product is XORed
into all its target slots at once and erased by X-basis measurement (Toffoli-free).
"""
from functools import lru_cache
import numpy as np
import kickmix as km
from generic_dqi import GenericDecoder

GF = km.GF2Field(12)


@lru_cache(None)
def cmat(c):
    """12x12 bit matrix of y -> c*y: column j is c * x^j."""
    return [[(GF.mul(c, 1 << j) >> i) & 1 for j in range(12)] for i in range(12)]


def gf_solve_vandermonde_inv(points):
    """Inverse of V[j][k] = a_j^k over GF(4096), by Gauss-Jordan."""
    n = len(points)
    M = [[GF.pow(a, k) if (a or k) else 1 for k in range(n)] + [1 if r == j else 0 for r in range(n)] for j, a in enumerate(points)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col])
        M[col], M[piv] = M[piv], M[col]
        inv = GF.invert(M[col][col])
        M[col] = [GF.mul(inv, v) for v in M[col]]
        for r in range(n):
            if r != col and M[r][col]:
                f = M[r][col]
                M[r] = [v ^ GF.mul(f, w) for v, w in zip(M[r], M[col])]
    return [row[n:] for row in M]  # Vinv[k][j]


@lru_cache(None)
def plan(n):
    """For deg-(n-1) y: points, eval coefs E[j][k] = a_j^k, and output coefs Out[p][k] (product p -> c_k)."""
    pts = [0] + [GF.pow(GF.primitive_element, j) for j in range(n - 1)]
    E = [[GF.pow(a, k) if (a or k) else 1 for k in range(n)] for a in pts]
    Vinv = gf_solve_vandermonde_inv(tuple(pts))
    out = []
    for j in range(n):
        out.append([Vinv[k][j] for k in range(n)] + [0])
    # p_inf contributes c_n directly and -sum_j Vinv[k][j] a_j^n to low coefficients.
    inf = []
    for k in range(n):
        s = 0
        for j, a in enumerate(pts):
            s ^= GF.mul(Vinv[k][j], GF.pow(a, n) if a else 0)
        inf.append(s)
    out.append(inf + [1])
    return pts, E, out


def lin_xor(b, pairs):
    """pairs: list of (c, src, dst): dst ^= c * src for registers."""
    srcs, dsts = [], []
    for c, src, dst in pairs:
        if not c:
            continue
        m = cmat(c)
        for i in range(12):
            for j in range(12):
                if m[i][j]:
                    srcs.append(src[j]); dsts.append(dst[i])
    if srcs:
        b.cx(km.array(srcs), km.array(dsts))


class ToomDecoder(GenericDecoder):
    def _poly_mul_xor(self, b, q1, q0, ys, targets):
        """targets[k] ^= (q0 + q1 x) * y, coefficient k (len(targets) == len(ys) + 1)."""
        n = len(ys)
        pts, E, out = plan(n)
        p = b.alloc_qubits(12)
        ya = b.alloc_qubits(12)
        qa = b.alloc_qubits(12)
        for j, a in enumerate(pts):
            ev = [(E[j][k], ys[k], ya) for k in range(n)]
            lin_xor(b, ev)
            lin_xor(b, [(1, q0, qa), (a, q1, qa)])
            b.ixor_gf2_mul(qa, ya, target=p, field=self.gf)
            lin_xor(b, [(out[j][k], p, targets[k]) for k in range(n + 1)])
            b.del_gf2_mul(qa, ya, target=p, field=self.gf)
            lin_xor(b, [(1, q0, qa), (a, q1, qa)])
            lin_xor(b, ev)
        b.ixor_gf2_mul(q1, ys[n - 1], target=p, field=self.gf)
        lin_xor(b, [(out[n][k], p, targets[k]) for k in range(n + 1)])
        b.del_gf2_mul(q1, ys[n - 1], target=p, field=self.gf)
        b.free(p); b.free(ya); b.free(qa)

    def _update(self, b, X, Y, i, q1, q0):
        T = 2 * self.t
        # Remainder: r_{i-1} = Y[0..T-i] -> X[0..T+1-i].
        self._poly_mul_xor(b, q1, q0, [Y[k] for k in range(T - i + 1)], [X[k] for k in range(T - i + 2)])
        # Cofactor: s_{i-1} degree d at Y[T-d] (d < i) -> degree d at X[T-d] (d <= i).
        self._poly_mul_xor(b, q1, q0, [Y[T - d] for d in range(i)], [X[T - d] for d in range(i + 1)])


def build_circuit(m=4095, n=70, chien=True):
    from fast_dqi import FastBuilder
    b = FastBuilder(km.CircuitBuilder())
    syn = [b.create_quantum_register(12, f'S{i}') for i in range(n)]
    mb = [b.create_classical_register(12, f'E{i}') for i in range(m)]
    ToomDecoder(n, km.GF2Field(12)).build(b, syn, mb, chien=chien)
    return b.finish_circuit()
