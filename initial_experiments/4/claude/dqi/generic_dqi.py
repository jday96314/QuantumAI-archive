"""Fixed-schedule ("generic") extended-Euclid Reed-Solomon MBUC decoder.

For random error patterns every quotient of the key-equation Euclid run on
(x^{2t}, S(x)) has degree exactly one, so step i (1..t) is
    r_i = r_{i-2} + (q1 x + q0) r_{i-1},   s_i = s_{i-2} + (q1 x + q0) s_{i-1}
with every coefficient at a classically known slot. Two 71-slot arrays hold
[r | s]: array A starts as (x^{2t}, 0), array B as (S, 1); odd steps update A,
even steps update B. The two remainder slots that step i zeroes become the
cofactor's two new top slots (cofactor degree d lives at slot 2t - d).

Step i runs only while `active[i]` (all previous steps active and the leading
coefficient of r_{i-1} is nonzero). With w errors this stops exactly after w
steps. A shot fails only when a leading coefficient vanishes by accident
(about w / 4096 per shot), which the repetition factor prices.

q is computed from the remainders before the update and erased from the
cofactors afterwards (q1 = s_i[i] / s_{i-1}[i-1], ...), so each step is reversible.
The result (sigma = s_w, omega = r_w) is in A when w is odd: one controlled swap
moves it to B, then the Chien search with Forney kickback runs on B.
"""
import kickmix as km
from fast_dqi import FastBuilder
from forney_dqi import product_terms


class GenericDecoder:
    def __init__(self, n_syndromes, field):
        self.n = n_syndromes
        self.t = n_syndromes // 2
        self.gf = field

    # ---- helpers -------------------------------------------------------
    def _gated_inverse(self, b, src, active):
        """Returns (g, inv, scaffold) with g = active * src^{-1}."""
        inv, scaffold = b.init_gf2_inverse_with_scaffold(src, field=self.gf)
        g = b.alloc_qubits(12)
        b.init_and(active, inv, g)
        return g, inv, scaffold

    def _del_gated_inverse(self, b, src, active, g, inv, scaffold):
        b.del_and(active, inv, g)
        b.free(g)
        b.del_gf2_inverse_with_scaffold(src, target=inv, scaffold=scaffold, field=self.gf)
        b.free(scaffold)
        b.free(inv)

    def _mul(self, b, x, y, target):
        b.ixor_gf2_mul(x, y, target=target, field=self.gf)

    def _q_from_r(self, b, X, Y, i, g, q1, q0):
        """q1 ^= X[top] g ;  q0 ^= (X[top-1] + q1 Y[lead-1]) g."""
        T = 2 * self.t
        self._mul(b, X[T + 1 - i], g, q1)
        tmp = b.alloc_qubits(12)
        b.cx(X[T - i], tmp)
        if T - i - 1 >= 0:
            self._mul(b, q1, Y[T - i - 1], tmp)
        self._mul(b, tmp, g, q0)
        if T - i - 1 >= 0:
            self._mul(b, q1, Y[T - i - 1], tmp)
        b.cx(X[T - i], tmp)
        b.free(tmp)

    def _q_from_s(self, b, X, Y, i, g2, q1, q0):
        """q0 ^= (s_i[i-1] + q1 s_{i-1}[i-2]) g2 ;  q1 ^= s_i[i] g2."""
        T = 2 * self.t
        tmp = b.alloc_qubits(12)
        b.cx(X[T - i + 1], tmp)
        if i >= 2:
            self._mul(b, q1, Y[T - i + 2], tmp)
        self._mul(b, tmp, g2, q0)
        if i >= 2:
            self._mul(b, q1, Y[T - i + 2], tmp)
        b.cx(X[T - i + 1], tmp)
        b.free(tmp)
        self._mul(b, X[T - i], g2, q1)

    def _update(self, b, X, Y, i, q1, q0):
        T = 2 * self.t
        for k in range(T + 1):
            if k <= T + 1 - i:  # remainder terms, Y remainder slots 0..T-i
                if 0 <= k - 1 <= T - i:
                    self._mul(b, q1, Y[k - 1], X[k])
                if k <= T - i:
                    self._mul(b, q0, Y[k], X[k])
            if k >= T - i:  # cofactor terms, degree d = T - k <= i
                d = T - k
                if 0 <= d - 1 <= i - 1:
                    self._mul(b, q1, Y[T - (d - 1)], X[k])
                if d <= i - 1:
                    self._mul(b, q0, Y[T - d], X[k])

    def _active_flag(self, b, Y, i, prev, flag):
        """flag ^= prev and (Y[lead] != 0)."""
        lead = Y[2 * self.t - i]
        z = b.alloc_qubits(1)[0]
        b.flip_if_equal(lead, 0, target=z)
        b.x(z)
        if prev is True:
            b.cx(z, flag)
        else:
            b.ccx(prev, z, flag)
        b.x(z)
        b.flip_if_equal(lead, 0, target=z)
        b.free(z)

    def _step(self, b, X, Y, i, active, forward):
        T = 2 * self.t
        q1 = b.alloc_qubits(12)
        q0 = b.alloc_qubits(12)
        r_lead = Y[T - i]
        s_lead = Y[T - i + 1]
        if forward:
            g, inv, sc = self._gated_inverse(b, r_lead, active)
            self._q_from_r(b, X, Y, i, g, q1, q0)
            self._del_gated_inverse(b, r_lead, active, g, inv, sc)
        else:
            g, inv, sc = self._gated_inverse(b, s_lead, active)
            # Inverse of _q_from_s: q1 first, then q0.
            self._mul(b, X[T - i], g, q1)
            tmp = b.alloc_qubits(12)
            b.cx(X[T - i + 1], tmp)
            if i >= 2:
                self._mul(b, q1, Y[T - i + 2], tmp)
            self._mul(b, tmp, g, q0)
            if i >= 2:
                self._mul(b, q1, Y[T - i + 2], tmp)
            b.cx(X[T - i + 1], tmp)
            b.free(tmp)
            self._del_gated_inverse(b, s_lead, active, g, inv, sc)
        self._update(b, X, Y, i, q1, q0)
        if forward:
            g, inv, sc = self._gated_inverse(b, s_lead, active)
            self._q_from_s(b, X, Y, i, g, q1, q0)
            self._del_gated_inverse(b, s_lead, active, g, inv, sc)
        else:
            g, inv, sc = self._gated_inverse(b, r_lead, active)
            # Inverse of _q_from_r: q0 first (needs q1), then q1.
            tmp = b.alloc_qubits(12)
            b.cx(X[T - i], tmp)
            if T - i - 1 >= 0:
                self._mul(b, q1, Y[T - i - 1], tmp)
            self._mul(b, tmp, g, q0)
            if T - i - 1 >= 0:
                self._mul(b, q1, Y[T - i - 1], tmp)
            b.cx(X[T - i], tmp)
            b.free(tmp)
            self._mul(b, X[T + 1 - i], g, q1)
            self._del_gated_inverse(b, r_lead, active, g, inv, sc)
        b.free(q1)
        b.free(q0)

    def _eea(self, b, A, Bv, flags, forward):
        t = self.t
        steps = range(1, t + 1) if forward else range(t, 0, -1)
        for i in steps:
            X, Y = (A, Bv) if i % 2 else (Bv, A)
            prev = True if i == 1 else flags[i - 2]
            if forward:
                self._active_flag(b, Y, i, prev, flags[i - 1])
            self._step(b, X, Y, i, flags[i - 1], forward)
            if not forward:
                self._active_flag(b, Y, i, prev, flags[i - 1])

    # ---- Chien search with Forney kickback ----------------------------
    def _chien(self, b, Bv, error_mbuc):
        gf = self.gf
        t = self.t
        T = 2 * t
        alpha = gf.primitive_element
        ainv = gf.invert(alpha)
        terms = product_terms(gf)
        sig = [Bv[T - d] for d in range(t + 1)]  # sigma_d
        om = [Bv[d] for d in range(t)]  # omega_d
        acc_s = b.alloc_qubits(12)
        acc_d = b.alloc_qubits(12)
        acc_o = b.alloc_qubits(12)
        inv = b.alloc_qubits(12)
        u = b.alloc_qubits(12)
        parity = b.alloc_qubits(1)[0]
        ladder = b.alloc_qubits(10)
        is_error = b.alloc_bits(1)[0]
        meas = b.alloc_bits(1)[0]

        def sums():
            for d in range(t + 1):
                b.cx(sig[d], acc_s)
                if d % 2:
                    b.cx(sig[d], acc_d)
            for d in range(t):
                b.cx(om[d], acc_o)

        def mask_map(mask):
            for k in range(12):
                b.push_condition(mask[k])
                for i_, j_ in terms[k]:
                    b.cx(inv[j_], u[i_])
                b.pop_condition()

        num = len(error_mbuc)
        for j in range(num):
            mask = error_mbuc[j]
            xj = gf.pow(ainv, j)  # evaluation point alpha^{-j}
            b.bit_store(is_error, True, control=mask)
            b.push_condition(is_error)
            sums()
            # e_j = omega(x) / sigma'(x) = x omega(x) / sum_{odd d} sigma_d x^d.
            if xj != 1:
                b.gf2_imul(xj, target=acc_o, field=gf)
            _, chain = b.init_gf2_inverse_with_scaffold(acc_d, target=inv, field=gf)
            mask_map(mask)
            for i_ in range(12):
                b.ccx(acc_o[i_], u[i_], parity)
            mask_map(mask)
            b.x(acc_s)
            b.init_and(acc_s[0], acc_s[1], ladder[0])
            for k in range(1, 11):
                b.init_and(ladder[k - 1], acc_s[k + 1], ladder[k]) if k < 10 else None
            b.ccz(ladder[9], acc_s[11], parity)
            for k in reversed(range(1, 10)):
                b.del_and(ladder[k - 1], acc_s[k + 1], ladder[k])
            b.del_and(acc_s[0], acc_s[1], ladder[0])
            b.x(acc_s)
            b.hmr(parity, meas)
            b.push_condition(meas)
            b.gf2_phase_by_product(acc_o, inv, mask=mask, field=gf)
            b.pop_condition()
            b.del_gf2_inverse_with_scaffold(acc_d, target=inv, scaffold=chain, field=gf)
            b.free(chain)
            if xj != 1:
                b.gf2_idiv(xj, target=acc_o, field=gf)
            sums()
            b.pop_condition()
            b.bit_store(is_error, False)
            # Advance coefficient registers: c_d *= alpha^{-d}.
            if j + 1 < num:
                for d in range(1, t + 1):
                    b.gf2_imul(gf.pow(ainv, d), target=sig[d], field=gf)
                for d in range(1, t):
                    b.gf2_imul(gf.pow(ainv, d), target=om[d], field=gf)
        if num > 1:
            for d in range(1, t + 1):
                b.gf2_idiv(gf.pow(ainv, d * (num - 1) % 4095), target=sig[d], field=gf)
            for d in range(1, t):
                b.gf2_idiv(gf.pow(ainv, d * (num - 1) % 4095), target=om[d], field=gf)
        b.bit_store(meas, False)
        b.free(km.array([is_error, meas]))
        b.free(acc_s + acc_d + acc_o + inv + u + km.array([parity]) + ladder)

    # ---- top level ----------------------------------------------------
    def build(self, b, syndromes, error_mbuc, chien=True):
        t = self.t
        T = 2 * t
        A = [b.alloc_qubits(12) for _ in range(T + 1)]
        Bv = [syndromes[k] for k in range(T)] + [b.alloc_qubits(12)]
        flags = b.alloc_qubits(t)
        b.x(A[T][0])
        b.x(Bv[T][0])
        self._eea(b, A, Bv, flags, forward=True)
        par = b.alloc_qubits(1)[0]
        b.cx(flags, par)
        for k in range(T + 1):
            b.cswap(par, A[k], Bv[k])
        if chien:
            self._chien(b, Bv, error_mbuc)
        for k in range(T + 1):
            b.cswap(par, A[k], Bv[k])
        b.cx(flags, par)
        b.free(par)
        self._eea(b, A, Bv, flags, forward=False)
        b.x(Bv[T][0])
        b.x(A[T][0])
        b.free(flags)
        b.free(Bv[T])
        for a in A:
            b.free(a)


def build_circuit(m=4095, n=70, chien=True):
    b = FastBuilder(km.CircuitBuilder())
    syn = [b.create_quantum_register(12, f'S{i}') for i in range(n)]
    mb = [b.create_classical_register(12, f'E{i}') for i in range(m)]
    GenericDecoder(n, km.GF2Field(12)).build(b, syn, mb, chien=chien)
    return b.finish_circuit()


if __name__ == '__main__':
    import sys
    c = build_circuit()
    print('qubits', c.num_qubits, 'toffolis', c.operation_statistics().num_toffoli_like)
    if len(sys.argv) > 1:
        c.to_file(sys.argv[1], format='kmb')
