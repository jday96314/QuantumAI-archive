# Follow-up research notes, not implemented gains

The current actual checkpoint is9.285435820940434. This file records avenues
considered after the angle-gauge improvement so they are not mistaken for
already-working circuits or repeated without checking earlier evidence.

The primary sources consulted for structured moment transforms were
[Sokal, A simple algorithm for expanding a power series as a continued fraction](https://arxiv.org/html/2206.15434),
[Comer and Kaltofen, Hankel matrices and Berlekamp/Massey](https://users.cs.duke.edu/~elk27/bibliography/10/CoKa10.pdf),
and generalized Schur/Toeplitz literature. Sokal's refined recurrence avoids
series inversion by retaining two adjacent series and taking scaled shifted
differences. That alone does **not** prove an N-field reversible implementation:
the adjacent series generally overlap in information but still occupy two
banks unless a further compressed reversible schedule is derived. Likewise,
Hankel symmetry does not by itself establish a2n-field generalized-Schur state;
the shift displacement can have a kernel requiring boundary data. No such
low-workspace moment-transform circuit has been built or resource-verified.
Do not substitute a speculative O(n²) product count for the tested compact
Euclid numbers inBLOCKED_EUCLID_RESEARCH.md.

Another conceptual route factors the locator into explicit field roots before
looking up runtime classical masks. A12-bit-address,12-bit-output arbitrary
runtime classical table can be represented in ANF, split6+6: generate114
nontrivial monomials, and use at most63*12=756bilinear products to coherently
load its output. The classical coefficient transform is linear. Its phase
unlookup can reuse the monomial banks. **This is an arithmetic proposal, not
an implemented or checked circuit.** Even if useful, an in-place reversible
locator factorization and magnitude evaluation are unresolved. Adding a fresh
420-qubit root bank can erase the Toffoli saving through the DQI Q² factor.

Do not repeat tree/projection probes as if they were new. Read-only experiment3
CONTINUATION.md records extensive Huffman/greedy/permutation tree searches
(up to30000proposals), failed8-bit500/600/700dictionary models, R2column
clustering, R4block fits, and backward pruning. Those outcomes constrain those
specific initializations, not all new tree or dictionary algorithms. The exact
signed gauge inANGLE_GAUGE_RESEARCH.md is distinct: it changes only encoding,
preserves every projector, and has full native verification.

The public Coulomb circuit inandigu/uncompute uses a different absolute-error
contract (32units at2^-24), not the competition's24-bit relative tolerance.
Its128-line resource count cannot be copied into this invsqrt task without a
new exhaustive accuracy and resource check. The existing local direct linear
screens already missed the selected quadratic circuit's cost/accuracy tradeoff.
