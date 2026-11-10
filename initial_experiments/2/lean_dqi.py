"""Narrower ToomDecoder.

* A 6-bit counter of active steps replaces the 35 per-step flags. Step i is
  active iff counter == i - 1 and lead(r_{i-1}) != 0. The flag is erased as
  (counter == i) after the controlled increment. The counter's low bit is the
  parity that selects which array holds the result.
* Evaluation points are formed in place: q(a) in q0 and y(a) in y_0.
* Gated inverses are deleted right after gating and recomputed to erase the gate.
* The Chien sums accumulate in place into sigma_0, sigma_1 and omega_0.
"""
import kickmix as km
from toom_dqi import ToomDecoder, plan, lin_xor
from forney_dqi import product_terms


class LeanDecoder(ToomDecoder):
    # ---- gated inverse without keeping the scaffold alive ---------------
    def _gate(self, b, src, active):
        inv, sc = b.init_gf2_inverse_with_scaffold(src, field=self.gf)
        g = b.alloc_qubits(12)
        b.init_and(active, inv, g)
        b.del_gf2_inverse_with_scaffold(src, target=inv, scaffold=sc, field=self.gf)
        b.free(sc); b.free(inv)
        return g

    def _ungate(self, b, src, active, g):
        inv, sc = b.init_gf2_inverse_with_scaffold(src, field=self.gf)
        b.del_and(active, inv, g)
        b.free(g)
        b.del_gf2_inverse_with_scaffold(src, target=inv, scaffold=sc, field=self.gf)
        b.free(sc); b.free(inv)

    def _r_relation(self, b, X, Y, i, g, q1, q0, compute):
        """compute: q1 ^= X[top] g, then q0 ^= (...) g.  erase: reverse order."""
        T = 2 * self.t
        def q0_part():
            tmp = b.alloc_qubits(12)
            b.cx(X[T - i], tmp)
            self._mul(b, q1, Y[T - i - 1], tmp)
            self._mul(b, tmp, g, q0)
            self._mul(b, q1, Y[T - i - 1], tmp)
            b.cx(X[T - i], tmp)
            b.free(tmp)
        if compute:
            self._mul(b, X[T + 1 - i], g, q1); q0_part()
        else:
            q0_part(); self._mul(b, X[T + 1 - i], g, q1)

    def _s_relation(self, b, X, Y, i, g, q1, q0, compute):
        """compute: q1 ^= s_i[i] g, then q0.  erase: q0 first, then q1."""
        T = 2 * self.t
        def q0_part():
            tmp = b.alloc_qubits(12)
            b.cx(X[T - i + 1], tmp)
            if i >= 2:
                self._mul(b, q1, Y[T - i + 2], tmp)
            self._mul(b, tmp, g, q0)
            if i >= 2:
                self._mul(b, q1, Y[T - i + 2], tmp)
            b.cx(X[T - i + 1], tmp)
            b.free(tmp)
        if compute:
            self._mul(b, X[T - i], g, q1); q0_part()
        else:
            q0_part(); self._mul(b, X[T - i], g, q1)

    def _step(self, b, X, Y, i, active, forward):
        T = 2 * self.t
        q1 = b.alloc_qubits(12)
        q0 = b.alloc_qubits(12)
        r_lead, s_lead = Y[T - i], Y[T - i + 1]
        first, second = (r_lead, self._r_relation), (s_lead, self._s_relation)
        if not forward:
            first, second = second, first
        g = self._gate(b, first[0], active)
        first[1](b, X, Y, i, g, q1, q0, True)
        self._ungate(b, first[0], active, g)
        self._update(b, X, Y, i, q1, q0)
        g = self._gate(b, second[0], active)
        second[1](b, X, Y, i, g, q1, q0, False)
        self._ungate(b, second[0], active, g)
        b.free(q1); b.free(q0)

    # ---- in-place evaluation points -------------------------------------
    def _poly_mul_xor(self, b, q1, q0, ys, targets):
        n = len(ys)
        pts, E, out = plan(n)
        p = b.alloc_qubits(12)
        for j, a in enumerate(pts):
            ev = [(E[j][k], ys[k], ys[0]) for k in range(1, n)]
            lin_xor(b, ev)
            lin_xor(b, [(a, q1, q0)])
            b.ixor_gf2_mul(q0, ys[0], target=p, field=self.gf)
            lin_xor(b, [(out[j][k], p, targets[k]) for k in range(n + 1)])
            b.del_gf2_mul(q0, ys[0], target=p, field=self.gf)
            lin_xor(b, [(a, q1, q0)])
            lin_xor(b, ev)
        b.ixor_gf2_mul(q1, ys[n - 1], target=p, field=self.gf)
        lin_xor(b, [(out[n][k], p, targets[k]) for k in range(n + 1)])
        b.del_gf2_mul(q1, ys[n - 1], target=p, field=self.gf)
        b.free(p)

    # ---- counter-based activity ------------------------------------------
    def _xor_activity(self, b, Y, i, c, a):
        """a ^= (c == i - 1) and lead(r_{i-1}) != 0."""
        e = b.alloc_qubits(1)[0]
        z = b.alloc_qubits(1)[0]
        b.flip_if_equal(c, i - 1, target=e)
        b.flip_if_equal(Y[2 * self.t - i], 0, target=z)
        b.x(z)
        b.ccx(e, z, a)
        b.x(z)
        b.flip_if_equal(Y[2 * self.t - i], 0, target=z)
        b.flip_if_equal(c, i - 1, target=e)
        b.free(z); b.free(e)

    def _eea(self, b, A, Bv, c, forward):
        t = self.t
        for i in (range(1, t + 1) if forward else range(t, 0, -1)):
            X, Y = (A, Bv) if i % 2 else (Bv, A)
            a = b.alloc_qubits(1)[0]
            if forward:
                self._xor_activity(b, Y, i, c, a)
                self._step(b, X, Y, i, a, True)
                b.iadd(1, target=c, control=a)
                b.flip_if_equal(c, i, target=a)
            else:
                b.flip_if_equal(c, i, target=a)
                b.isub(1, target=c, control=a)
                self._step(b, X, Y, i, a, False)
                self._xor_activity(b, Y, i, c, a)
            b.free(a)

    # ---- Chien with in-place sums -----------------------------------------
    def _chien(self, b, Bv, error_mbuc):
        gf = self.gf
        t = self.t
        T = 2 * t
        ainv = gf.invert(gf.primitive_element)
        terms = product_terms(gf)
        sig = [Bv[T - d] for d in range(t + 1)]
        om = [Bv[d] for d in range(t)]
        acc_d, acc_s, acc_o = sig[1], sig[0], om[0]
        odd = [sig[d] for d in range(3, t + 1, 2)]
        even = [sig[d] for d in range(2, t + 1, 2)]
        src = km.array([q for r in odd + [sig[1]] + even + om[1:] for q in r])
        dst = km.array([q for r in [acc_d] * len(odd) + [acc_s] * (1 + len(even)) + [acc_o] * (t - 1) for q in r])
        rsrc = km.array(list(reversed(list(src))))
        rdst = km.array(list(reversed(list(dst))))
        inv = b.alloc_qubits(12)
        u = b.alloc_qubits(12)
        parity = b.alloc_qubits(1)[0]
        ladder = b.alloc_qubits(10)
        is_error = b.alloc_bits(1)[0]
        meas = b.alloc_bits(1)[0]

        def mask_map(mask):
            for k in range(12):
                b.push_condition(mask[k])
                b.cx(km.array([inv[j_] for i_, j_ in terms[k]]), km.array([u[i_] for i_, j_ in terms[k]]))
                b.pop_condition()

        num = len(error_mbuc)
        for j in range(num):
            mask = error_mbuc[j]
            xj = gf.pow(ainv, j)
            b.bit_store(is_error, True, control=mask)
            b.push_condition(is_error)
            b.cx(src, dst)
            if xj != 1:
                b.gf2_imul(xj, target=acc_o, field=gf)
            _, chain = b.init_gf2_inverse_with_scaffold(acc_d, target=inv, field=gf)
            mask_map(mask)
            for i_ in range(12):
                b.ccx(acc_o[i_], u[i_], parity)
            mask_map(mask)
            b.x(acc_s)
            b.init_and(acc_s[0], acc_s[1], ladder[0])
            for k in range(1, 10):
                b.init_and(ladder[k - 1], acc_s[k + 1], ladder[k])
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
            b.cx(rsrc, rdst)
            b.pop_condition()
            b.bit_store(is_error, False)
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
        b.free(inv + u + km.array([parity]) + ladder)

    def build(self, b, syndromes, error_mbuc, chien=True):
        t = self.t
        T = 2 * t
        A = [b.alloc_qubits(12) for _ in range(T + 1)]
        Bv = [syndromes[k] for k in range(T)] + [b.alloc_qubits(12)]
        c = b.alloc_qubits(6)
        b.x(A[T][0]); b.x(Bv[T][0])
        self._eea(b, A, Bv, c, forward=True)
        for k in range(T + 1):
            b.cswap(c[0], A[k], Bv[k])
        if chien:
            self._chien(b, Bv, error_mbuc)
        for k in range(T + 1):
            b.cswap(c[0], A[k], Bv[k])
        self._eea(b, A, Bv, c, forward=False)
        b.x(Bv[T][0]); b.x(A[T][0])
        b.free(c); b.free(Bv[T])
        for a_ in A:
            b.free(a_)


def build_circuit(m=4095, n=70, chien=True):
    from fast_dqi import FastBuilder
    b = FastBuilder(km.CircuitBuilder())
    syn = [b.create_quantum_register(12, f'S{i}') for i in range(n)]
    mb = [b.create_classical_register(12, f'E{i}') for i in range(m)]
    LeanDecoder(n, km.GF2Field(12)).build(b, syn, mb, chien=chien)
    return b.finish_circuit()
