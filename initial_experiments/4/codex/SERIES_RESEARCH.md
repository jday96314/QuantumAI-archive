# Low-workspace polynomial arithmetic

These are research circuits. None is part of the verified submission.

`series_inverse.py` implements normalized formal power-series inversion modulo
z^n in characteristic two, with three scratch qubits. The constant term must
be one. It composes power maps with exponents 3, 5, 17, ... . Their cumulative
exponents are 2^(2^r)-1, which become -1 in the truncated unit group.

For a power map F -> F^(1+q), split F=A+z^m B with m=ceil(n/2).
The high part becomes B*A^q + [A^(1+q)]_high; the low part recursively becomes
A^(1+q). Truncated unit triangular products are reversible in place. Middle
products use native-field Toom evaluation and transposition, with linear
Frobenius maps. Targets and sources are disjoint in every field product.

| Length | Q | Toffolis | Quantum check |
| --- | ---: | ---: | --- |
| 35 | 423 | 12,804 | 512 shots, PASS |
| 70 | 843 | 37,494 | 512 shots, PASS |

The optional `--optimized` variant selects direct sparse products when fewer
than the Toom points and bounds short coefficient arrays. Its separate
artifacts/checks pass512shots at length35 **423Q/11166T**, and length70
**843Q/32622T**. The original results above remain preserved. The new Möbius
series map and compact Euclid blocks use this arithmetic; see
BLOCKED_EUCLID_RESEARCH.md. These numbers do not replace a whole-decoder count.

`check_series_inverse.py` compares every output coefficient against an
independent convolutional reciprocal, checks the constant term, and checks
all scratch and phase. It includes zero-tail and all-one series. Evidence:
`results/series_inverse_35_check.json`, `results/series_inverse_70_check.json`.

This primitive makes a serial continued-fraction decoder possible within
approximately the syndrome storage, but its complete cost is not attractive.
The sum of inversion costs at lengths 70,68,...,2 is **517,818 Toffolis** in one
direction. Reversing this and adding the current scan already exceeds 1.38
million gates, before normalizing coefficients or reconstructing locator and
numerator. The inherited complete decoder uses 645,755. Lower width may repay
some of this, but there is no demonstrated score improvement. Repeated full
reciprocals should not be integrated without a substantially faster schedule.

`build_bezout_quantum.py` is a separate structured linear-map experiment.
For generic degree-17 E,O, a reversible Euclidean factorization applies
Bezout(E,O) to a 17-word payload, restores E,O, and uses 46 extra qubits and
49,624 Toffolis. The source Euclid stores its low quotient directly in g[16],
avoiding a 12-qubit temporary. A full Sylvester map for E[0]=1,
H=E*V+O*z*U, uses 46 extra qubits and 62,098 Toffolis. Both quantum checks pass
512 generic random cases. Degree-gap pairs are explicitly rejected by these
tests, so neither is a complete decoder. The two Sylvester directions alone
cost 124,196 gates, slightly exceeding the 122,850-gate theoretical scan saving
from replacing an 80-gate kernel with a 50-gate kernel at 4,095 points.

Source attribution for the reused field arithmetic and Euclidean generator:
Apache-2.0 `andigu/uncompute`, revision
64b5a13aa249f1360846657392e6c66817b26d0f, retained read-only under
`../../3/vendor/uncompute`. Native-field prototypes must be converted through
the established isomorphism before use with host GF(4096) syndrome inputs.
