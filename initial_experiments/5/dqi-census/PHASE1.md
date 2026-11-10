# Decoder architecture census: phase 1, before component construction

Status: analytical hypotheses, not verified decoders. The incumbent remains selected.
All estimates below charge forward, phase extraction, and reversal; no quantum
branch is discounted by its occurrence probability. Q is maximum physical width.

Read: round-5 repair, phase, hybrid, core, dispatch RESULTS; round-4 REPORT and
Codex continuation; Khattar et al. arXiv:2510.10967v1; public discovered README.
The public README is a different cost objective (TQ), with a fixed phase mask.
It does not certify a gain under this host's R Q² T and per-shot masks. The
leaderboard's 18.9 concerns the product of three task ratios, not a DQI certificate.

Reference: Q961, mean T1019865.856484375, R1.6968460732519504,
C1598204229384.9265. Core649802, scan370152 static. Zero failures in a full
102400-shot run would imply R1.1912604648, not R1. Classical correctness is
not evidence of that native outcome. Define P=4024.017333984375 active scan
points after the independently checked 70-point interpolation (70 guaranteed,
additional random zero masks in expectation). The old compact unpadding
lifetime U=9792 is charged unless a replacement is explicitly supplied.

Every candidate has two R cases: R0=1.1912604648 CONDITIONAL on zero native
failures, and Rref=1.6968460733 if the incumbent exception behavior survives.
The R0 objectives in the ranking are research projections, never accepted scores.
The present compact core actually has Rref; using R0 requires a new charged
gap-correct core. Any such additional T and Q is an explicit unknown, G.

## A. Direct reciprocal-product end state, linear scan

Decode to the existing compact (B,A,pinv,completion history) in the 70 syndrome
words. After unpadding, monic B has degree L and A is the scaled prior cofactor.
Replace A IN PLACE by F'=(A B')^-1 mod B. At B(beta)=0 the magnitude is
beta^70 F'(beta); absorb beta^70 into each classical trace mask, so evaluation
of B,F' and the root test cost 11T/point. The standard evaluator's different
normalization must not be substituted here. Reverse the same in-place map,
repad, and reverse decoding. On the split squarefree, unit-A domain the map
is an involution: f(f(A))=A. This is an algebraic identity, not a free erasure.

Peak table (target layout, not allocated):

| Register | Q |
|---|---:|
| B,A,scale in pinned syndrome storage | 840 |
| completion/history, peak | 35 |
| materialized leading word | 12 |
| temporary budget allowed at Q961 | 74 |
| Total target | 961 |

T=649802+U+11P+2C+G = 703858.1907+2C+G. C is the COMPLETE per-direction
map including inverse preparation, scratch erasure, and exceptional inputs.
At C=100000,G=0,Q961,R0: about 0.994e12 (1.61x). This is an allowance,
not an algorithmic estimate. Even C=0 only projects 2.06x at Q961: this
architecture also needs a cheaper/narrower core to reach 2.5x.
It beats the same-core 88T scan iff C<154924.6674 when both scans use interpolation; without
interpolation the familiar threshold is 157657.5 per direction. Both include
forward AND reverse. Against the actual un-interpolated incumbent mean scan, the threshold is158003.9154. The U term cancels in these comparisons.
Failures: existing degree gaps plus nonunit/repeated-root or bad completion
handling in the converter. Domain restriction must be reversible; Rref until
repaired. Biggest uncertainty: can the reciprocal-product involution use <=74
extra Q and <~157kT/direction? Test the quotient-algebra Frobenius route first,
with actual quantum-modulus reduction, then distinguish that route from EEA.

## B. Lazily allocated BM windows, then compact storage

Process an early prefix using BM's coefficient windows. Allocate only the
statically touched sigma/template slots; quantum degree never controls physical
allocation. Consumed syndrome words host the BM template. At round r=2k,
convert the processed prefix and active template to a compact cofactor/remainder
state, erase the external locator and BM history, and continue compact updates.
A handoff must be an invertible coordinate change; copying the compact state
and running BM backward is a separately priced wide alternative. End with
(B,A), use the 88T scan, reverse compact suffix, inverse handoff, reverse BM.
A is optional further work, not silently included.

Peak table at r=20 from existing exact lazy-BM replay, before handoff:

| Register | Q |
|---|---:|
| Pinned syndrome/template storage | 840 |
| Materialized locator and remaining persistent private state | 251 |
| temporary pool (maximum over prefix) | 54 |
| Prefix peak | 1145 |
| Any handoff scratch beyond reusable pool | H_Q (unknown) |

The existing prefixes already allocate lazily. This proposal adds a genuine
representation transition; it does not claim the already-used laziness as a
new saving. Let B20 be exact forward+reverse prefix T, E10 the compact
forward+reverse suffix from degree10, and H the one-way transition T:
T=B20+E10+2H+370152; Q>=1145 before new scratch. A planning value T=900000
would cost about1.406e12 at R0, only1.14x; exact prefix/suffix attribution must
replace it. Later handoff reduces suffix T but raises Q (r30:1200,r60:1380).
Failures: BM windows/stack, handoff zero diagonals, existing compact gaps;
Rref is a comparison scenario and R0 conditional; neither bounds new handoff errors. Biggest uncertainty: existence of
an in-place prefix-to-compact map, including history erasure. A full-width
late BM core cannot be substituted while claiming the early-prefix peak.

## C. Syndrome-resident locator AND evaluator throughout (rational coordinates)

Represent processed data by sigma(0)=1, omega=sigma S mod x^r, and unprocessed
syndrome tail. On a valid degree-k processed prefix, omega[k:r]=0. Put the
locator in those established zero rows, rather than storing another 432Q
locator. Maintain a rational correction descriptor in the remaining processed
slots, and allow a small retained defect/parking record when a leading
coefficient vanishes. At completion the 70 words hold sigma and omega (35+35
variable coefficients; fixed constant handled by wiring), so the 80T Forney
scan consumes exactly the state the decoder produces. Reverse the rational
updates to recover the original syndrome. No reciprocal-derivative bank.

Peak target table:

| Register | Q |
|---|---:|
| Resident sigma/omega/active tail/correction descriptors | 840 |
| Defect/parking/history allowance | 41 |
| Scalar/update scratch allowance | 80 |
| Total target | 961 |

At r=2k there are k zero rows, enough for k variable locator coefficients;
this dimension count does NOT prove the correction descriptor also fits.
A nonminimal reference update divides V=sigma S by old sigma, changes sigma
by alpha*x^m*tau, and multiplies V by new sigma. It exposes the aliasing issue:
storing sigma inside V makes ordinary disjoint-register triangular arithmetic
invalid. Native component must measure this reference and test whether
protected zero rows can stay untouched. It must not clear an input by promise
when its location contains sigma.

Planning core K=280000 (near BM arithmetic227247 plus ~53k for resident
updates/control), scan80P, adapters10000: T=611921.3867, Q961, R0 -> ~0.673e12
(2.37x). K=240000,Q940 gives ~0.602e12 (2.66x). These are explicit targets,
not sums of implemented components. At Q961 the 2.5x budget permits only
259162T core+adapter after the scan. Biggest uncertainty: a cheap alias-safe
reversible rational update with its correction/history in 70 words. Failures:
coefficient support, singular corrections, and lifetime errors; no intrinsic
RS decoding failures for distinct nonzero errors <=35. R0 only after all are
handled; otherwise a native-measured R must replace it. Test the reference
coordinate update, including exact work and the co-live register obstruction.

## D. Structured Hankel / PGZ via displacement elimination

Build H_ij=S_(i+j+1) implicitly, solve the Prony recurrence by rank-revealing
structured elimination, produce sigma, then omega via a unit triangular map.
A rank-deficient H is NORMAL for weights <35, so fixed full-rank inversion
cannot be the decoder. Pivot/rank history is retained for reversal. Scan80P,
restore syndromes by inverse elimination. This is linear-system elimination,
not an EEA/BM schedule.

An honest general Cauchy-like displacement layout has two length-35 generators
on EACH side (140 field words), not two vectors total. Fixed left/right nodes
are classical. A symmetric characteristic-two representation needs additional
diagonal information; symmetry does not magically remove it.

| Register | Q |
|---|---:|
| Original syndrome / generator data, share 70 words if overwritten | 840 |
| Other two displacement generators | 840 |
| output locator, if not yet folded into dead generator rows | 420 |
| pivots/rank/scratch planning allowance | 100 |
| Peak general layout | 2200 |

Optimistic arithmetic 3*35² field products/direction ->308700T lifetime,
scalar inversions2*35*66=4620, control allowance20000, scan80P: T655241.3867.
At Q2200,R0 cost3.778e12. A novel symmetric compressed generator layout
might reduce Q to~1100 ->0.944e12, but missing diagonals and reversible
rank changes are not free. Dense Gaussian PGZ is much worse (35² field
matrix=14700Q before other state; cubic T). Biggest uncertainty is a
rank-revealing compressed generator layout, NOT whether ordinary Gaussian
elimination works. Retain as a research family; the explicit wide layout loses.

## E. Wide BM + evaluator + Psi

Use Gu's locator recurrence with its stack; monicize, overwrite syndromes with
omega and 35 zero tail words; tail becomes the derivative-inverse pass; scan
23P; reverse Psi, product, normalization, BM. This spends width for a large T
cut and is fully priced under Q².

| Register | Q |
|---|---:|
| Syndromes/omega/tail workspace | 840 |
| Locator | 432 |
| BM template | 444 |
| Euclid replay pool | 48 |
| BM reversal stack | 19 |
| extra cofactor word + one word | 24 |
| discrepancy + degree counter | 20 |
| Euclid event record + preswap | 7 |
| emitter pivots | 3 |
| Published whole-circuit peak | 1837 |

Published T442397 is a fixed-mask count. Runtime masks: 4024.0173 mean active
points rather than4024; add~0.399T to its point term before host conversion.
At R0,Q1837,T442397.4 ->1.778e12, a 0.90x ratio. R is not public-paper R:
cutoff/stack/Euclid failures must be natively measured under host confidence.
The known syndrome-consuming BM variation reduces width but its actual Q1429
uses an80T scan; giving it a23T scan adds a new tail-lifetime problem.
Biggest uncertainty: lifetime-compatible removal of ~500Q at equal T, not
MBU's already-established reverse saving. The pinned upstream liveness intervals give this co-live partition at first
1837Q peak, Toffoli175139, in euclid/R; see results/published_peak.json.

## F. Root-list / amplitude decoder with indexed phase

Find sigma, factor its split squarefree polynomial to a coherent list of <=35
roots, and recover amplitudes. Store (root,value) pairs. Apply each input
classical mask by a quantum-addressed lookup, then reverse amplitudes,
factorization, and syndrome decoding. This changes phase extraction: no Chien
sweep. Root order must be canonical or its permutation history retained.

| Register | Q |
|---|---:|
| Pinned original syndromes | 840 |
| 35 root indices and 35 magnitudes | 840 |
| sigma if retained during factoring | 420 |
| factor/rank/order scratch allowance | 100 |
| Peak straightforward layout | 2200 |

A unary-iteration 4096-address lookup costs ~4094T PER ROOT (classical data
does not make a quantum address free): 143290T for35 lookups before lookup
cleanup, factorization, or decoder. Even granting full decode/factor/value
lifetime300000T, T>=443290 and R0,Q2200 ->2.556e12. Compact in-place roots
would require a new syndrome representation: at Q1100 this optimistic lower
projection is0.639e12, barely2.5x before omitted costs. Biggest uncertainty:
in-place, canonical factorization with no retained 840Q syndrome copy.
Failures include irreducible/non-squarefree intermediates, ordering and rank;
R0 requires all promised split locators and arbitrary masks. No free random
factoring success assumption is allowed.

## G. Transposed residue / transform-domain phase

Use the residue identity sum_(B(beta)=0) M(beta)/[A(beta) B'(beta)]
= coefficient_(x^(L-1)) of M/A mod B, with M incorporating beta^70 and the
per-shot trace-dual mask interpolant. Classically interpolate all4095 masks.
Quantumly reduce/contract by the rational denominator, apply one trace phase,
then reverse the reduction and decoder. This avoids explicit error values.

| Register | Q |
|---|---:|
| Compact B,A and scale | 840 |
| residue/reduction accumulator,35 words | 420 |
| degree/history plus scalar scratch allowance | 100 |
| Streaming peak | 1360 |

A simple streaming reduction needs about4095*35=143325 field products
per direction: T>=2*42*143325+649802=12689102, even ignoring quotient
retention. Cost atQ1360,R0 >27e12. Fast transposed remainder trees can
trade this time for many live coefficient banks. Dynamic classical M makes
multiplication by M linear, but reciprocals of quantum B are not classical.
Biggest uncertainty: a small-memory transposed modular reduction beating
~100kT/direction. This route is distinct from point scan preparation.
Failures inherit the compact gaps and singular residue denominator; Rref until
fixed. The elementary streaming realization loses; fast transposed contraction
remains open, not disproved.

## H. Wronskian state / inverse-only scan

Replace one compact cofactor by W=A B'+A' B. At roots W=A B', so the scan
uses the inverse-only46T cubic-trace body, no42T point product. Decode directly
into (B,W), or reversibly convert the compact pair; reverse that coordinate
map after the scan. The characteristic-two W is a square polynomial, so it
has only~35 independent field coefficients when degree<=68.

| Register | Q |
|---|---:|
| Compact input/state B,W and scale | 840 |
| completion/history | 35 |
| adapter workspace allowance | 86 |
| Target (unallocated) | 961 |

Measured prior squarefree repadding7140T/direction; generic Sylvester map
62098T/direction (does not handle gaps). Planning conversion lifetime
2*(7140+62098)=138476, plus core649802+46P ->973382.7974T.
AtQ961,R0 cost1.071e12 (1.49x), before any missing H-division/total gap
adapter. A direct W state could remove conversion costs, but would need a
new update law. Biggest uncertainty: total, low-workspace W state transition
and H-division, not the already-tested padding. Rref until gaps are solved.

## Phase-1 ordering and decisions

Conditional objective ordering of concrete targets: C~0.673e12, A~0.994e12,
H~1.071e12, B~1.406e12, E~1.778e12, F>=2.556e12, D~3.778e12,
G>27e12. C and A are selected for minimal native probes. C's expected cost
has the largest uncertainty and no implementation: its ranking is a
resource-target ranking, not a defensible posterior mean. A's100k conversion
allowance likewise is not measured. No scalar expected value can honestly
be assigned before resolving those uncertainties. Updated measurements must
replace these placeholders, and a failed schedule must not kill the family.

No candidate is dropped here just because its first estimate loses. Explicit
losing realizations are recorded. Unresolved architectures remain open. The
next phase will test C's rational-coordinate update with physical liveness,
and A's modular Frobenius cost including quantum-modulus reduction. A lazy
BM prefix attribution will also sharpen B. If a resulting design projects
>=2x, a concrete gate/lifetime/failure/verification build plan follows.
