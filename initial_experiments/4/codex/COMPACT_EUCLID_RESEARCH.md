# Compressed Euclidean decoder state

These components are excluded from the submission. They do not constitute a
complete decoder and do not improve its verified score.

Later work implements quantum generic steps, quotient blocks, full armed
sequences and completion padding. See **BLOCKED_EUCLID_RESEARCH.md** for their
tests, exact resources, current cost limitation and the next H-map adapter.

`compact_euclid.py` retains two growing cofactors and the high band of one
remainder. For x^N and S, the invariant `P B + R Bprev = x^N` determines omitted
low remainder coefficients. After k generic degree-one quotient steps, retain
`R[k:N-k]`. This state occupies N+1 field words. Its forward and backward maps
match ordinary Euclid through every tested stage and round-trip exactly.
`results/compact_euclid_check.json` records the cases and skipped degree gaps.

`monic_compact_euclid.py` normalizes both cofactors to monic form and retains
one inverse leading-remainder scalar. The invariant is

    P B + R Bprev = p x^N.

Implicit monic leading coefficients leave exactly N stored field words. With
r=lead(R), s=sublead(R), and pinv=1/p, the next step uses

    beta = r pinv
    v = s/r
    alpha = v + sublead(B)
    Bnew = (x+alpha) B + beta Bprev
    newpinv = 1/r.

The new high remainder band follows from `B Rnew = r x^N + Bnew R`. The test
checks every remainder and cofactor against full Euclid, every single-step
inverse, and full round trips. It passes 370 generic cases spanning N4..70;
degree gaps remain outside this diagnostic's scope. The serial arithmetic
estimate is **411,516 field-product Toffolis per direction**, before inversions,
cleanup, controls, gap handling, or the phase scan. This estimate does not justify
building a whole serial decoder; a cheaper blocked update is still needed.

`monic_euclid_first.py` implements the first step in native quantum gates. For
a=S[N-1] and v=S[N-2]/a, it keeps the high remainder values
`S[j-1]+v S[j]`, v, and 1/a in the input array. A retained flag arms a zero
pivot reversibly. At N70 it passes 512 shots for forward values, inverse values,
one-builder composition, phase, and scratch cleanup. It uses **875 qubits** and
**3012 Toffolis per direction**. N8 also passes at 131 qubits and 408 Toffolis.
See `results/monic_first_70_check.json` and `results/monic_first_8_check.json`.

`inverse_mask_kernel.phase_inverse` supplies the ungated inverse-parity phase
needed to erase an inverse's old source. It uses **36 Toffolis**, with 34 total
qubits including the 12 quantum input bits. The standalone check covers all
4096 field elements against each of twelve mask basis vectors, plus zero and
random masks: **57,344 shots PASS**, restoring both public operands and zeroing
scratch. See `results/inverse_mask_ungated_check.json`.

The classical compressed models use the attributed upstream field convention;
the quantum first-step and mask primitive use the host's native field convention.
Each test compares values within its own field. Source attribution and the
Apache-2.0 upstream revision are documented in README.md.
