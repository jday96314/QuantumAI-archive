"""LeanDecoder with a lower qubit peak.

* q is allocated after the first gate, and erased by X-basis measurement
  before the second gate is formed. The erasure phase is applied by
  Toffoli-free gf2_phase_by_product fixups; only the trilinear part of q0
  needs one recomputed product (q1 Y) and its deletion.
* The Chien zero-test ladder reuses u's qubits after the parity is formed.
"""
import kickmix as km
from lean_dqi import LeanDecoder
from forney_dqi import product_terms


class Lean2Decoder(LeanDecoder):
    def _erase_by_measurement(self, b, g, q1, q0, top, nxt, ynext):
        """Erase q1 = top g and q0 = (nxt + q1 ynext) g given gate g."""
        m1 = b.alloc_bits(12); m0 = b.alloc_bits(12)
        b.hmr(q1, m1); b.hmr(q0, m0)
        b.free(q1); b.free(q0)
        b.gf2_phase_by_product(top, g, mask=m1, field=self.gf)
        b.gf2_phase_by_product(nxt, g, mask=m0, field=self.gf)
        if ynext is not None:
            r1 = b.alloc_qubits(12)
            self._mul(b, top, g, r1)
            tmp = b.alloc_qubits(12)
            self._mul(b, r1, ynext, tmp)
            b.gf2_phase_by_product(tmp, g, mask=m0, field=self.gf)
            b.del_gf2_mul(r1, ynext, target=tmp, field=self.gf)
            b.free(tmp)
            b.del_gf2_mul(top, g, target=r1, field=self.gf)
            b.free(r1)
        b.bit_store(m1, False); b.bit_store(m0, False)
        b.free(m1); b.free(m0)

    def _step(self, b, X, Y, i, active, forward):
        T = 2 * self.t
        r_lead, s_lead = Y[T - i], Y[T - i + 1]
        if forward:
            first, second = (r_lead, self._r_relation), s_lead
            erase = (X[T - i], X[T - i + 1], Y[T - i + 2] if i >= 2 else None)
        else:
            first, second = (s_lead, self._s_relation), r_lead
            erase = (X[T + 1 - i], X[T - i], Y[T - i - 1])
        g = self._gate(b, first[0], active)
        q1 = b.alloc_qubits(12); q0 = b.alloc_qubits(12)
        first[1](b, X, Y, i, g, q1, q0, True)
        self._update(b, X, Y, i, q1, q0)
        # Erase q: measure, then the phase fixups need the second relation's gate.
        m1 = b.alloc_bits(12); m0 = b.alloc_bits(12)
        b.hmr(q1, m1); b.hmr(q0, m0)
        b.free(q1); b.free(q0)
        self._ungate(b, first[0], active, g)
        g = self._gate(b, second, active)
        top, nxt, ynext = erase
        b.gf2_phase_by_product(top, g, mask=m1, field=self.gf)
        b.gf2_phase_by_product(nxt, g, mask=m0, field=self.gf)
        if ynext is not None:
            r1 = b.alloc_qubits(12)
            self._mul(b, top, g, r1)
            tmp = b.alloc_qubits(12)
            self._mul(b, r1, ynext, tmp)
            b.gf2_phase_by_product(tmp, g, mask=m0, field=self.gf)
            b.del_gf2_mul(r1, ynext, target=tmp, field=self.gf)
            b.free(tmp)
            b.del_gf2_mul(top, g, target=r1, field=self.gf)
            b.free(r1)
        self._ungate(b, second, active, g)
        b.bit_store(m1, False); b.bit_store(m0, False)
        b.free(m1); b.free(m0)


def build_circuit(m=4095, n=70, chien=True):
    from fast_dqi import FastBuilder
    b = FastBuilder(km.CircuitBuilder())
    syn = [b.create_quantum_register(12, f'S{i}') for i in range(n)]
    mb = [b.create_classical_register(12, f'E{i}') for i in range(m)]
    Lean2Decoder(n, km.GF2Field(12)).build(b, syn, mb, chien=chien)
    return b.finish_circuit()
