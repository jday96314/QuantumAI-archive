# Phase 1 — architecture census, before building

2026-10-24. Incumbent 182956.5271606425; 1.5× target 121971.0181;
1.3× planning cutoff 140735.7901. These are hypotheses, not built resources.
T below means executed Toffolis; synthesis is additional. No selection changes.

Read-only evidence: inv-width, inv-phase, inv-dyadic RESULTS and source. In
particular, removing a named bank does not necessarily remove physical wires:
the incumbent borrows phase headroom for product bits and carries. All peak
budgets below include replacement arithmetic scratch, and need compilation.

For a numerical margin m turns and R arbitrary rotations use
S(R,m) = R*(0.56*log2(R/(pi*m))+5.3)/4 + direct_T/4.
The unchanged incumbent polynomial has m=19995185/2^64 = 1.0839411508124785e-12.
S(34,m)=250.85021 including its one direct T. New fits are budgeted at
m=1e-12 unless otherwise specified. A fit does not establish that margin.

## A — conditionally clean input as output headroom (top prototype)

A normalized input M has its lowest k bits zero, while the output expands by
s=floor(k/2). Borrow ten input bits, arranged in reverse order, as the high
bits of a 37-position cyclic phase bank. Its last s positions are guaranteed
zero. Rotation by s therefore inserts zeros into the low end of the polynomial.
Only the first s head positions belong to the polynomial: phase high position
j only if s>j. Phase the low 27 positions unconditionally. Restore the rotation
before restoring the input. This differs from dirty-input arithmetic borrowing:
no padding subtraction is needed, and borrowing occurs only during phasing.
Allocate arithmetic scratch on demand, rather than reserving headroom globally.

| Peak ownership | Q budget |
|---|---:|
| normalized input (also borrowed phase headroom) | 22 |
| exponent flags | 5 |
| affine residual | 23 |
| updated slope | 17 |
| retained curvature | 0–4 |
| carry counters, row products, adder workspace | 9–19 |
| total | 76–90 |

T hypothesis 1750–2150, R about 37 (formerly Clifford high phases become
controlled non-Clifford phases), m unchanged; S about 274. Cost 154k–218k.
An aggressive release/replay schedule at Q=70,T=1730 would cost 140.3k and
qualify for a build plan, but is **not** in the conservative range above.
Biggest uncertainty: can the arithmetic fit below 80 without consuming the
savings, and what is the real cost of occupancy-controlled phasing? Exhaust
borrowed-word semantics; then compile resource schedules with the new lifetime.

## B — exponent through a multiplexed rotation (second prototype)

Keep only the unshifted 27-bit polynomial. For each bit, apply its phase with
angle selected by s=0..10. An exact Walsh expansion of the exponent selector
uses parity-controlled Z rotations, with CNOT parity formation. Extend the
function on unreachable s=11..15 to minimize the rotation count. This removes
both the output shift and its permanent headroom. It is distinct from phasing
each Booth row and paying a separate shift for each row.

| Peak ownership | Q budget |
|---|---:|
| input | 22 |
| exponent | 5 |
| constant residual | 23 |
| slope | 17 |
| retained curvature | 0–4 |
| temporary products/carries and parity wire | 9–19 |
| total | 76–90 |

Optimistic T=1565 (217 shift gates removed with no arithmetic penalty).
If R<=35, S about259 and cost at Q=76 is138.6k (1.32x); if R=200, S about1550
and cost exceeds236k. m unchanged. Biggest uncertainty: sparse *exact* phase
synthesis on the eleven reachable exponent values. Count exact Walsh support,
including controlled-phase decomposition, before claiming this projection.

## C — self-delimiting in-place exponent

Pack k into zeros in M. For k=0 retain M. For k>=1 use
E=M-2^21+2^(k-1). Then k=ctz(E)+1, and M is recovered by replacing the marker
with the fixed high one. This is injective and uses 22 rather than 27 persistent
wires. Decode parity/scale controls only when needed. The local-coordinate
marker must be removed or absorbed into the fit; treating E as M is incorrect.

| Peak ownership | Q budget |
|---|---:|
| packed input including exponent | 22 |
| transient decoded controls | 1–3 |
| constant/slope/curvature | 23+17+4 |
| headroom and arithmetic workspace | 19–22 |
| total | 86–91 |

T=1850–2250, R=34, target m=1e-12, S about251; objective181k–228k.
Combination with A could reach aboutQ=72, T=1800, cost148k. Biggest uncertainty:
repeated prefix decode and marker repair, which can erase a five-qubit saving.
Packing is an information representation result, not a free arithmetic operand.

## D — exponent-scaled coefficients with a fixed physical binary point

Fit each unnormalized dyadic input interval directly to phase units, so exponent
scaling is in the coefficient values. A fixed approximately36-bit phase word
replaces the 27+10 shifting bank; output shift gates disappear. Coefficient
banks can be emitted in slices and erased by measurement with phase repair.

| Peak ownership | Q budget |
|---|---:|
| original input (no normalized copy) | 22 |
| interval/coordinate routing | 4 |
| physical phase accumulator | 36 |
| slope | 18 |
| curvature | 7 |
| products and carries | 12 |
| total | 99 |

T=1700–2600 with exponent-subtree sharing assumed; R=33; m target1e-12,
S about244; objective192k–282k. Biggest uncertainty: table sharing between
exponents; an explicit elevenfold table is prohibitive. This eliminates a
shift but does not make its ten bits disappear from numerical precision.

## E — cubic on coarser aligned intervals

Use a cubic on 32 intervals per parity instead of the current 128–256 adaptive
quadratics. Fit a rounded/truncated pipeline jointly, with nominal fractional
precisions 27/19/11/4. Release each high coefficient as soon as the next Horner
stage consumes it, replaying only for reverse arithmetic. Affine-precondition c0.

| Peak ownership | Q budget |
|---|---:|
| input and flags | 27 |
| affine c0 and final slope | 23+18 |
| surviving higher coefficients | 0–4 |
| headroom/products/carries | 22 |
| total | 90–94 |

T=1750–2350 (table about62, but three products each direction), R=34,
target m=1e-12, S about251; objective180k–245k. Biggest uncertainty: whether
coefficient rounding permits the tiny cubic coefficient and short inner products.
Fewer table leaves alone do not justify an improvement claim.

## F — Gu-style dyadic linear with compressed coefficient loading

Power-of-two-aligned leaves, a cyclic shared local coordinate, in-place
slope-to-destination accumulation, and measured table erase with phase repair.
At relative2^-24, expect roughly1900–2048 parity-augmented line entries rather
than the paper's128 for an absolute2e-6 target. New twist: split high coefficient
bits and small low residuals into separate QROM slices, or a two-row SelectSwap
lookup; recycle the selected blocks before arithmetic.

| Peak ownership | Q budget |
|---|---:|
| input and scale | 27 |
| constant and slope | 27+18 |
| one coefficient slice / routing | 8–16 |
| arithmetic workspace (slice storage recycled) | 10 |
| total | 90–98 |

T=2100–3000 even assuming roughlyhalf ordinary QROM load, R=34,
target m=1e-12, S about251; objective212k–319k. Biggest uncertainty: sharing
lookup work without retaining a second full coefficient word. inv-width's
uncompressed linear failure is evidence about that lookup, not this compression.

## G — exact rational quadratic with fused square/middle product

The existing coefficients evaluated as
c0 - c1*z/2^12 + c2*z^2/2^23 have exact all-input margin1.6008221224300678e-11.
Fuse truncated square columns directly into a narrow destination; do not form
a full square or whole c2*z^2 product. Erase selected monomials with measured
AND repair. This differs from the prior direct circuit's exact wide square.

| Peak ownership | Q budget |
|---|---:|
| input and exponent | 27 |
| one coefficient bank (phased/released separately) | 23 |
| narrow destination / residual | 18 |
| streamed square/product/carries | 12 |
| total | 80 |

Ambitious T=1300–1800, R=50–90, S about349–640, cost132k–195k.
At T=1400,R=50 it projects140k and qualifies for a build plan. Biggest uncertainty:
column truncation and coefficient replay can cost more than the saved product;
50 rotations assumes shared exponent routing. Neither bound is built evidence.

## H — two-iteration table-seeded Newton correction

Use a short seed and exploit a small residual 1-M*y0^2. Store the residual at
reduced precision and directly add -y0*residual/2 to the phased seed. The seed
must be refitted to the actual omitted columns; no full general-width Newton.

| Peak ownership | Q budget |
|---|---:|
| input and exponent | 27 |
| seed | 13 |
| small residual | 14 |
| final correction | 24 |
| product/carry workspace | 14 |
| total | 92 |

T=2000–3200, R=55–70, target m=1e-12, S about415–534;
objective222k–344k. Biggest uncertainty: producing a narrow residual without
holding a full y0^2 or losing relative accuracy. Not the earlier serial Horner.

## Ranking and next phase

A and B test a common bottleneck without changing the numerical approximation,
so their margin is known and failures can be attributed to the new mechanism.
B has a conditional1.3× projection; A is the wider, lower-rotation alternative.
G is the strongest unresolved arithmetic architecture, conditional on a new
fused primitive; its Phase 3 plan must state those assumptions explicitly.
No architecture is accepted or globally ruled out by these estimates.

Gu reference: https://arxiv.org/html/2610.22367v1, Supplement III.1–III.3.
The techniques are adapted, not the paper's looser error budget or gate totals.
No upstream Apache-2.0 code has been copied at this stage.
