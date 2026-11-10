Full-circuit follow-up: the trimmed decoder is now ACCEPTED; see the acceptance block at the top of RESULTS.md. The research-only acceptance status below is historical.

# Compact Euclid: cost attribution and limited reduction

The best prototype reduces the core from **649,802 to 624,206 T (3.939%)**
at the selected **961-Q** peak. It misses the requested 20% reduction and the
852-Q target. The width-constrained schedule search returns exactly the
selected block boundaries. Under the requested <10%-reducible stopping clause,
this investigation stops without building or accepting a replacement decoder.
This is evidence about this compact representation and the tested arithmetic,
not a lower bound on other decoders or a claim that the leaderboard is saturated.

## Reference and evidence

The reference is round-4 Codex
`candidates/cubic_trace_product_n70_m4095.kmb`, SHA-256
`4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0` in its
existing native report. That report records 102,400 shots, zero value errors,
380 phase errors, Q=961, mean T=1,019,865.856484375, host-confidence R=
1.6968460732519504 and DQI objective 1,598,204,229,384.9265. These are inherited
verification results, not new verification of this prototype. Reference files
were never modified. Overall reference 9.3677221531 is preserved.

`attribute_core.py` replays the selected core's actual first transform,
block3/suffix/monomial completion, serial-after-22, virtual zero, lazy flags,
low-workspace inversion and final cofactor scaling. Product-term and direct
builder gate accounting reproduces **324,901 forward + 324,901 reverse**.
The method suppresses Clifford emission, so a whole prototype KMB is not implied.
The 12 native block circuits independently match the corresponding counted
components exactly; `summarize.py` asserts these comparisons.

## Attribution

These categories are exclusive and sum to 324,901 in **each** direction.
Scaling and inverse gadgets inside cofactor/head operations appear in their
own rows; they must not be added again from the inclusive totals below.

| Operation / arithmetic | Forward T | Reverse T |
| --- | ---: | ---: |
| Series multiplication/division and middle products | 162,582 | 162,582 |
| Cofactor updates and serial boundary products | 58,074 | 58,074 |
| Field scaling | 53,466 | 53,466 |
| Controlled monomial padding | 26,244 | 26,244 |
| Inverses, including first/serial replacement chains | 16,320 | 16,320 |
| Six-moment Möbius transfer, remaining products | 1,176 | 1,176 |
| Six-moment continued fraction, remaining products | 2,940 | 2,940 |
| Flag/completion equality tests | 385 | 385 |
| First remainder shear and remaining completion preparation | 3,714 | 3,714 |
| **Total** | **324,901** | **324,901** |

The head Möbius maps cost **23,184 T** per direction inclusive of their
scaling, inverses and triangular products. Head continued-fraction transforms
cost **13,335 T** inclusive. Together, six-moment transfer is **36,519 T**,
not an additional cost on top of the table. Serial cofactor regions also
contain tail transforms; their inclusive 150,480 T is likewise not additive.

The middle-product attribution further separates the actual operations:

| Region | Forward T | Reverse T |
| --- | ---: | ---: |
| Tail prefix-preserving suffix multiplication / its inverse | 60,564 | 60,564 |
| Tail ordinary series division / corresponding reverse multiply | 64,302 | 64,302 |
| Tail ordinary series multiply / corresponding reverse divide | 20,076 | 20,076 |
| Head triangular products and divisions | 17,640 | 17,640 |
| **Total** | **162,582** | **162,582** |

The selected core's 50.04% series share is the main cost. Its convolution
emitter uses evaluated/transposed bilinear products rather than emitting every
schoolbook pair, so counting dense coefficient pairs overstates possible gains.

## Measurement erasure and reverse asymmetry

| Forward temporary | Present handling | Additional reverse saving |
| --- | --- | --- |
| Product fan-out bit, bilinear evaluation scratch | Measured and phase-repaid within each product | None by a second erasure |
| Inversion N/M/S/Y chain | Early erasure; selected low-workspace helper recomputes phase-repair data | Selective retention saves 36 T per eligible inverse, in both directions |
| Monic scaling product word | Product replaces the old word; old word is measurement-erased | None without changing the state transform |
| Six-word copied block head | Copy/transform/erase completed inside each block | No head remains at the scan/reverse boundary |
| Boundary old-sublead and completion preparation words | Locally computed, consumed and measurement-erased | No retained redundant history |
| 35 completion flags | Persist through the scan and govern unpadding | Their original pivot tests are no longer available without reconstructing prior states |
| Two cofactors and inverse scale | The entire surviving 840-bit state | These encode consumed syndrome information and must be inverted to restore input |

The scan boundary retains only the compact 840-bit state and 35 completion
flags; all core-private scratch is clean. Thus there is **no additional
forward-history bank** that can be measurement-erased to make the existing
reverse selectively cheaper. The estimated *mechanical asymmetric-reverse*
saving for this architecture is **0 T**. This does not rule out changing its
representation or reverse algorithm.

Applying the supplied Gu ratio mechanically would give
324,901×1.61=523,090.61 T, a 19.5% core saving. That number is a counterfactual,
not a resource estimate for an implemented reverse. The current forward already
uses measurement erasure; equal directional totals do not mean that it uses
unitary uncomputation for every temporary.

## Prototypes

`components/trimmed_middle.py` trims coefficient-independent zero rows and
columns before choosing schoolbook versus evaluated middle products. For

    high[i] ^= sum_j A^q[L + offset + i - j] low[j],

if the highest nonzero coefficient index of A^q is d, only the final
min(L,d−offset) low words and the first min(H,d−offset) high words can matter.
After trimming, the same formula holds with the reduced L. Evaluation uses
the transposed coefficient positions shifted by `offset`; this also handles
partial suffix maps. All operand/target linear maps are Clifford operations,
and no additional quantum storage is allocated. The improvement saves
**11,214 T per direction**, 3.4515% of the core.

`components/budget_inverse.py` keeps the tower inversion chain when the actual
live count plus 30 inverse bits and 12 body-scratch bits fits within 961 Q.
The source remains unchanged during the body, so its chain can be erased for
free afterward. Otherwise it calls the selected low-workspace inverse. Exactly
**44 inverses per direction** qualify, saving **1,584 T per direction**.
Keeping all 146 inverse contexts would save only 5,256 T per direction and
would exceed the selected peak in the busy blocks. Even the deliberately
optimistic combination of trimming and making *all* 16,320 inverse T free is
only an 8.47% core reduction; the remaining products cannot be removed by
erasing inverse temporaries.

Combined: **312,103 forward + 312,103 reverse = 624,206 T**. Series work falls
to 151,368 T and inverse work to 14,736 T per direction; all other exclusive
categories remain unchanged. Native block widths include 961 Q at k=1 and
k=19, 920 Q at k=22, and 932 Q at k=34. Scan storage remains the selected
959-Q layout by unchanged interface/allocation analysis; the full circuit has
not been emitted or verified.

Native checks comprise **28 circuits × 512 shots = 14,336 component shots**:
16 trimmed middle-map circuits with q=1,2,4, selected-sized tails, nontrivial
suffix offsets, zero coefficients and both composition orders; and 12 N70
block/serial circuits at k=1,19,22,34, with independent polynomial outputs,
forward/inverse/composition, forced zero-tail completion, clean scratch and
zero phase. These sampled checks establish components, not full decoding.

## Algorithmic alternatives and stopping evidence

**Running remainder products / ratio recurrence.** The existing compact map
uses a numerator middle product followed by old-cofactor series division.
Maintaining both products independently would need an additional tail-sized
bank early. Replacing the tail by its full local ratio avoids the copied head
and three cofactor passes, but Schur recurrence now inverts and scales the
*whole surviving series*. `ratio_probe.py` tests that implementation with the
new trimmed products and retained inverses:

| Three-step entry | Full-ratio forward T | Q | Best local block3 forward T |
| --- | ---: | ---: | ---: |
| k=1 | 175,953 | 889 | 15,687 |
| k=10 | 145,758 | 898 | 32,769 |
| k=19 | 81,660 | 907 | 35,007 |
| k=31 | 15,291 | 919 | Serial remains cheaper in the schedule DP |

These are exact emitter counts, not additional native simulation results.
The experiment rejects this particular full-ratio recurrence, not every
running-product formulation. A representation that propagates omitted low
remainder coefficients cheaply remains open algebraic work.

**Larger early blocks / half-gcd.** Four-step blocks with the original schedule
reduce forward T to 307,837 (5.25%) but peak at 984 Q. Even a standalone early
four-step block with the new arithmetic peaks at 968 Q. Its copied eight-word
head exceeds the 961-Q allowance before late flags accumulate. This does not
rule out a different in-place half-gcd algorithm. `schedule_screen.py` counts
120 feasible serial/b1/b2/b3 local alternatives in both directions, then uses
dynamic programming over degree boundaries under Q≤961. The optimum remains
block3 at k=1,4,7,10,13,16,19, followed by serial k=22..34, for **624,206 T**.
No additional schedule gain was found.

**End state (sigma, omega) or (sigma, F).** The present final coordinates are
monic Bprev (34 explicit words), monic B (35), and pinv (1). The identity

    P B + R Bprev = p x^70

gives R = p x^70 / Bprev mod B on the coprime generic branch. Therefore an
in-place modular reciprocal could exchange (Bprev,pinv) for a 35-word
evaluator without enlarging the field bank. A total reversible circuit with
completion/gap handling at the requested width was not obtained. The current
35 retained flags alone give 875 live bits at the scan, already above 852;
changing endpoint coordinates alone does not remove them.

An 80-T Forney body versus the selected 88-T body saves only
8×4,095=**32,760 T** over the scan. At unchanged Q/R, an endpoint adapter would
need to cost less than **16,380 T per direction**, including all completion
handling, to pay for itself. A full decoder that naturally ends in these
coordinates is a separate design opportunity. No (sigma,omega) end state or
852-Q decoder is claimed here.

## Objective projection and acceptance

The unchanged scan adds 370,152 T, giving a projected static total of
**994,358 T**, versus 1,019,954 selected: **2.5095%** complete-circuit reduction.
Using static T as a mean-T approximation and the requested R≈1.2:

| Assumption | R Q² T |
| --- | ---: |
| Selected static T, Q=961, assumed R=1.2 | 1,130,338,725,640.8 |
| Prototype static T, Q=961, assumed R=1.2 | **1,101,972,593,421.6** |
| Prototype static T, hypothetical Q=852, assumed R=1.2 | 866,170,139,558.4 |

Q=852 is not achieved. R=1.2 is not measured; the selected measured R is 1.69685.
The actual acceptance objective remains host Clopper–Pearson
R(epsilon_plus)×Q²×mean_Toffolis. No full-default 102,400-shot run or independent
salts were performed because the counted gain is below the requested scope
threshold. No replacement is accepted, and no score improvement is claimed.

## Attribution and reproduction

The replay and liveness interpretation reuse round-5 `dqi-hybrid` and
`dqi-phase/baseline_audit.py`; product-term counting follows round-4
`profile_compact_parts.py`. Arithmetic is imported read-only from
andigu/uncompute, Apache-2.0, revision
64b5a13aa249f1360846657392e6c66817b26d0f. The original license remains in
`initial_experiments/3/vendor/uncompute/LICENSE`. New files and artifacts are
confined to this assigned directory. Reproduction commands are in RESULTS.md.
