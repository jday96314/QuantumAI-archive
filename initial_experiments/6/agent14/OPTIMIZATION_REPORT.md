# Quantum circuit optimization report

## Result and scope

The final complete submission achieved an **audited local score of 14.064304039685362**, exceeding the target of 14. The main advance was a different DQI decoder representation, followed by reductions in arithmetic and temporary quantum storage. Inverse square root improved through exact circuit changes that preserved the existing approximation. The chemistry experiments produced useful negative results, but no new chemistry component beat the inherited selection.

The comparison below uses the best separately verified components available from rounds 4 and 5. Their calculated combined score is 9.76167122781691; this comparison does not assume that those components had already been packaged together. Relative to that component baseline, the final combined score improved by about 44.1%.

| Problem | Starting task cost | Final task cost | Cost reduction in this run |
|---|---:|---:|---:|
| DQI | 1,435,930,717,449.94 | 515,978,772,816.41 | **64.07%** |
| Inverse square root | 182,956.53 | 170,242.25 | **6.95%** |
| Chemistry | 696,942,846.91 | 696,942,846.91 | **0%** |

Lower task costs are better. The overall score is the geometric mean of the three host-baseline-to-submission cost ratios. These task-specific costs therefore cannot be added together.

The final artifacts are [submission.zip](submission.zip) and [submission.summary.json](submission.summary.json). The archive SHA256 is `dbcc81d3dcd80e416f203ff63f3539db65347cdfd27636259f68192224ff3396`. This report describes completed results as of October 25, 2026; the chronological [research log](RESULTS.md) also contains older intermediate states.

## DQI

### Objective and starting point

DQI must recover error information from 70 quantum syndrome words over GF(4096), apply the required phase using 4,095 classical masks, and restore the syndrome register and scratch space. The input alone occupies 840 qubits. The score charges

`repetition factor × peak qubits² × mean Toffolis`.

This makes width especially important. Failures also matter: a smaller circuit can lose if its error rate requires more repetitions. The inherited decoder used 1,004 qubits and about 1,195,802 mean Toffolis, with zero observed failures in the starting default report.

### Changes that worked

**Replace the decoder representation.** The largest gain came from quadratic-extension Schur extraction and continued-fraction synthesis. Interpolation nodes lie outside GF(4096), in a quadratic extension. Their fixed node polynomials have the form `q_i(x) = x² + x + κ_i`.

The early implementation maintained affine frames, scalar normalization, and coherent node swaps. These components worked, but frame updates and exceptional-case bookkeeping were expensive. The later raw continued-fraction representation eliminated scalar normalization and node swapping from extraction. Its heads are `l_i(x) = a_i + b_i x`, with recurrence

`R_i = l_i + q_i / R_(i+1)`.

Synthesis then updates a numerator and denominator by

`(N, B) ← (l_i N + q_i B, N)`.

This representation reduced both arithmetic and stored bookkeeping. Degree-one polynomial products were evaluated and interpolated together, avoiding two independent field multiplications per input coefficient. The raw architecture produced the first complete submission above 13.

**Prepare an amplitude polynomial before scanning.** Instead of performing a field division at every scanned point, the circuit transforms the numerator into `C = A / B′ mod B`. At locator roots, the remaining known field factors can be absorbed into the classical masks. The scan then consists of polynomial evaluation, a zero test on the locator, and a masked phase on the amplitude.

This reduced the scan from 80 Toffolis per point in an earlier Schur implementation to **11 per point**, or 45,045 across all 4,095 points. The conversion has a substantial preparation and reversal cost, but it is still cheaper than repeating the expensive arithmetic at every point. Holding selected X-basis measurement outcomes across the scan also avoided keeping their quantum intermediates alive.

**Reduce the endpoint to three fields and choose a better scale.** The terminal continued-fraction state originally used four infinity-jet fields. A three-field state `(p, q, d)` carries the required information. For a head `(a, b)`, its update is

`p′ = q`, `q′ = p + bq`, `d′ = d + pq + (a + b)q²`.

For nonzero `q`, choosing the terminal tail `N = [d/q, p]`, `B = [q]` makes the final synthesized denominator monic. The physical zero-`q` case has a separate reversible representation. This choice removed both whole-polynomial normalization passes and simplified the outer fold. The complete circuit fell from **557,645 to 547,147 Toffolis** at the same 896-qubit peak.

**Consume a redundant remainder during quotient recovery.** The last large improvement changed reverse synthesis. At a quadratic node, the remainder obeys `num = head × den`. Once the quotient is recovered, retaining a separate copy of `num` is unnecessary.

The new circuit measures the remainder in the X basis and reuses its 24 wires for norm and inverse work. It repays the measurement phase using the bilinear identity above. Polynomial restoration proceeds in the evaluation basis, so the old remainder does not have to be recreated. The proper outer fold needs an explicit correction for its two absent high coefficients; including that correction was necessary for correctness.

The quotient kernel fell from **664 to 348 Toffolis**, with private workspace falling from 20 to 16 qubits. Across the complete circuit, this saved **12,572 Toffolis and four qubits**, producing the selected **892-qubit, 534,575-Toffoli** circuit.

**Improve finite-field primitives and their lifetimes.** Supporting changes included an emitter that forms operand parities inside existing registers, delayed repayment of measurement phases, and a compact inverse-parity kernel. The selected scalar inverse machinery combines 66-Toffoli creation with a 30-Toffoli phase cleanup. These changes were tested with values, scratch restoration, and phase all checked.

A separate 28-AND field-cubing primitive also beat generic 42-AND multiplication and improved an early accepted decoder. It was superseded by the final architecture. The final circuit continues to reuse the existing 42-AND general field multiplier; this run did not find a cheaper replacement.

### Final resource breakdown

These are static Toffoli counts from the selected component builds. They sum exactly to the complete circuit count; quantum workspace is reused between stages.

| Stage | Toffolis |
|---|---:|
| Forward raw extraction | 114,240 |
| Forward rational synthesis | 34,750 |
| Forward amplitude conversion | 80,412 |
| Phase scan | 45,045 |
| Reverse amplitude conversion | 82,662 |
| Reverse rational synthesis | 45,376 |
| Reverse raw extraction | 132,090 |
| **Total** | **534,575** |

The final peak is only 52 qubits above the 840-qubit input. Compared with the starting decoder, mean Toffolis fell by about 55.3%, and peak width fell from 1,004 to 892.

### Approaches that did not win

| Approach | Observed result | What it established |
|---|---|---|
| Early Schur decoder without pivot repair | Q949 and 885,278 Toffolis, but 756 value errors and one phase error; cost about 1.641 trillion | The lower arithmetic count did not offset the repetition penalty. |
| Coherent node swaps and affine frames | Complete verified implementations, but costly controls, normalization, and stored state | A valid intermediate architecture; the raw representation removed much of its overhead. |
| Blocked continued fractions | Generic algebra tests passed, and larger blocks reduced the measured tail cost | Tail-only counts omitted extraction, transfer matrices, cleanup, and exceptional cases; they did not demonstrate a better complete circuit. |
| Removing completion bookkeeping too aggressively | A prototype confused a genuine polynomial factor with an artificial completion marker | Correctness required information from the preceding polynomial. The corrected version was useful but later superseded. |
| Further field-multiplier searches | Different constructions and local searches retained 42 ANDs | No improvement in the searched families; this is not a lower-bound proof. |

The selected decoder still has rare unsupported extension-pole cases. Both the default and fixed independent full runs observed **five value errors and zero phase errors in 102,400 shots**. The verifier includes these errors in its confidence and repetition penalty; the result is not a claim of exact decoding on every input. The resource savings outweighed that penalty. Eight saved failure cases also behaved identically before and after the final synthesis changes under two fixed measurement seeds.

Evidence: [final build](results/raw_cf_destructive_n70_m4095_build.json), [default verification](results/raw_cf_destructive_n70_verify_default.json), [independent verification](results/raw_cf_destructive_n70_verify_independent.json), and [remainder-consuming quotient](destructive_quotient.py).

## Inverse square root

### Objective and inherited implementation

This problem requires the phase proportional to `1/√x` for every integer from 1 through 3,133,452, with the input restored and scratch clean. The score charges peak qubits times mean Toffolis plus the equivalent cost of rotation synthesis. A numerical approximation that barely meets the tolerance can therefore be expensive because little margin remains for synthesizing its rotations.

The inherited implementation already used binary normalization, a fitted quadratic approximation, affine coefficient preconditioning, radix-four Booth multiplication, and measured cleanup. It used **90 qubits**, **1,793 maximum Toffolis**, and 34 non-Clifford rotations. Its task cost was 182,956.53.

The selected improvements in this run preserve that approximation, its rounding behavior, and its phase margin.

### Changes that worked

**Share adjacent Booth products.** Two adjacent source bits multiplied by the two Booth magnitude controls form a short binary polynomial product. Four elementary products can be replaced by three:

`A = c₀s₀`, `B = c₁s₁`, `C = (c₀ + c₁)(s₀ + s₁)`.

The middle coefficient is `A + B + C` over GF(2). The existing Booth controls are disjoint, which permits the product combination used by this integer-arithmetic implementation. The optimization also preserves the original rounding and carry rules. Extending the pairing to the rounded carry boundary saved further gates.

The first complete paired version used Q89/T1727; the carry-aware version used **Q89/T1709** and passed all 200,540,928 prescribed shots. Its cost was **173,447.64**.

**Free the curvature register and reconstruct it later.** Keeping the curvature coefficient alive during the final linear multiplication contributed to the width peak. An exact quadratic Boolean decoder can reconstruct it from values already live at that point: 21 low constant-coefficient bits, eight coarse address bits, and four high slope bits.

The selected decoder uses 52 shared quadratic products plus an exception correction. This lets the curvature register be measured and freed, then recreated during reverse evaluation to repay its measurement phase. Workspace scheduling brings the complete circuit to **85 qubits and 1,763 maximum Toffolis**.

That is more Toffolis than the Q89/T1709 intermediate, but its lower width wins under the actual objective. The final cost is **170,242.25**, a **6.95% reduction** from the inherited circuit.

| Complete implementation | Qubits | Maximum Toffolis | Task cost |
|---|---:|---:|---:|
| Inherited quadratic implementation | 90 | 1,793 | 182,956.53 |
| Paired Booth products | 89 | 1,727 | 175,049.68 |
| Carry-aware product pairing | 89 | 1,709 | 173,447.64 |
| Curvature reconstruction | **85** | **1,763** | **170,242.25** |

All four rows have full prescribed-shot verification. The final circuit retains 34 rotations and the same worst error-to-tolerance ratio, about 0.97967.

### Approaches that did not win

| Approach | Observed result | Why it lost or failed |
|---|---|---|
| Narrower curvature coefficient | A numerical fit passed all inputs, but the best tested circuit was Q88/T1863, cost about 185,533 | More lookup leaves and arithmetic outweighed the width saving. |
| Fused rounded-square correction | Useful coordinate and square kernels passed exhaustive tests; the complete prototype reached Q90/T2703 | Separate coefficient lookup and replay, square preparation, arithmetic, and cleanup made the full implementation too expensive. |
| Coarser-address cubic Horner approximation | A fitted cubic passed numerical checks with 167 leaves; complete prototype Q92/T1982 | The smaller lookup did not pay for another polynomial stage and its live registers. |
| Radix-eight multiplication prototype | The initial implementation failed native checks | Magnitude-three partial products require integer carries; XORing them directly was invalid. |
| Coefficient bit embedding | Two curvature bits could be embedded with 324 leaves; tested three- and four-bit constraints failed accuracy | The successful case enlarged the lookup too much. |
| First curvature-reconstruction schedule | The standalone decoder passed, but the complete preflight found dirty scratch | It kept coefficient bits live while treating them as zero workspace. The selected decoder uses the actual live values and retains the original cleanup timing. |

The main lesson was to optimize the exact arithmetic and register lifetimes together. A cheaper isolated polynomial or squarer did not imply a cheaper complete phase circuit.

Evidence: [final full verification](results/pair_curve_live_q85_t1763_verify_default.json), [paired arithmetic](booth_pair_carry_adder.py), [curvature reconstruction](curve_live_restore.py), and [fused-prototype stage costs](results/fused_stage_profiles.json).

## Chemistry

### Objective and selected component

The chemistry score is `λ_eff × total Toffolis × total qubits`. A candidate must also satisfy independent one-body and two-body error limits of 20 mHa, the gap and normalization conditions, and the per-oracle qubit cap. Physical vectors must be evaluated after decoding their integer Givens-angle words.

The selected component was inherited from round 5: **one block, 719 dictionary vectors, 350 factor columns, eight-bit angles**, with structured four-word ANGLE batches. It uses **454 scored qubits**, **31,608 Toffolis**, and `λ_eff = 48.56733747`, for cost **696,942,846.91**. Its two-body error is about **19.989973 mHa**.

The earlier successful compiler composition reduced scored width from 478 to 454 while increasing Toffolis from 30,504 to 31,608, giving about a 1.58% cost reduction at unchanged normalization. That improvement belongs to the inherited round-5 work. This run retained it unchanged.

### Experiments and outcomes

**Repair smaller dictionaries directly on the integer lattice.** The most developed path started from a 671-vector subdictionary with error 23.943463 mHa. Greedy integer repair, full-angle Hessian proposals, and lattice proposals improved it, but did not cross the required 20 mHa boundary.

| B671 method | Best resulting two-body error |
|---|---:|
| Starting exact-word dictionary | 23.943463 mHa |
| Integer descent | 21.127105 mHa |
| Full-angle Hessian proposals | 21.105792 mHa |
| Lattice proposals | **21.105777 mHa** |

The final normalization was also worse, around 56.28. A global QUBO search, 15 rounds of eight-atom dictionary exchanges, and 88 nearby shared-tree changes did not improve that endpoint.

For those final fixed words, allowing the full 671-column positive-semidefinite coefficient fit still left **20.034993 mHa** error. Increasing coefficient rank alone therefore could not rescue that dictionary. This is a numerical capacity result for those words, not a limit on every possible B671 dictionary.

**Reduce factor rank while keeping the larger dictionary.** At B719, rank C340 could pass the tensor-error limit, with error about **19.910161 mHa**. However, the tested factor gauges increased normalization substantially. Forty L1 Jacobi sweeps reduced `λ_eff` only to **56.1249**, still well above 48.5673. The small preparation saving could not compensate.

**Use continuous fitting and quantization-aware methods.** Some smaller continuous geometries passed the tensor-error limit. For example, R2/B383/C200 reached about **19.379 mHa** before quantization. Rounding and refitting destroyed that margin; the tested ADMM continuation ended around **156 mHa**. A separate R4/B255/C90 continuation remained around **69 mHa**. Straight-through and projected straight-through optimization on B671 did not improve its starting integer checkpoint and often worsened it substantially.

These experiments exposed a representation problem: a good continuous fit can depend on directions that are poorly preserved by the eight-bit Givens grid. A passing floating-point residual was insufficient evidence for a usable native factorization.

**Split or overlap dictionaries.** Six disjoint partitions of the inherited dictionary failed even an unrestricted block-covariance span screen, with errors of roughly **281–1,274 mHa**. Twelve overlapping dictionaries also failed, at roughly **189–868 mHa**. These relaxed fits were useful rejection tests for the particular dictionaries, before paying for positivity, rank, or oracle construction.

Giving every physical factor its own small dictionary also failed with the tested spectral initializer: truncating to B7 left about **4.77 Ha** error. Keeping all 76 directions reproduced the floating-point incumbent, which checked the implementation but offered no useful small-block solution.

### What the negative results mean

Three requirements repeatedly conflicted: accuracy after integer decoding, low normalization, and a smaller resource product. Improving only one often worsened the others.

The inherited fixed-shift normalization bound is `λ_eff > 46.78618`. Compared with 48.56734, this leaves less than 4% normalization headroom even in principle under that model. A large chemistry improvement would therefore need lower circuit resources as well as a physically valid factorization. None of the tested alternatives achieved that combination, so no chemistry replacement was selected.

Evidence: [selected native chemistry report](../../5/chem-census/results/structured_native.json), [B671 rank capacity](results/rank_capacity_671.json), [rank-gauge experiment](results/rank_gauge_C340.json), [dictionary exchanges](results/exchange_B671_final.json), [tree search](results/tree_B671_final.json), and the inherited [normalization-bound derivation](../../4/codex/NORMALIZATION_BOUND.md).

## Verification and saved milestones

Selection required complete circuits and the unchanged native scoring rules. The final inverse circuit passed **200,540,928 shots**. DQI passed **102,400 default shots** and a separate **102,400-shot independent sample** with a declared fixed salt. Chemistry reused its unchanged **50,880-shot** native evidence, bound to the exact canonical chemistry digest.

The submission ZIP was reopened, extracted, and hashed. Its members were checked against the verified component digests. Earlier milestone archives were preserved, and all five milestone ZIP hashes were checked at completion.

| Milestone | First complete audited score |
|---|---:|
| 10 | 10.268368324127126 |
| 11 | 11.02830129721277 |
| 12 | 12.078134489004718 |
| 13 | 13.202469964520382 |
| 14 | **14.064304039685362** |

The most transferable result was that representation and cleanup dominated many local arithmetic gains. DQI improved most when normalization and redundant state disappeared from the algorithm. Inverse square root improved when existing information could reconstruct a coefficient and release its storage. Chemistry remained constrained by the joint behavior of integer geometry, normalization, and circuit resources.

The final audit records are [archive_audit.json](results/archive_audit.json) and [all_milestones_final_audit.json](results/all_milestones_final_audit.json). The first archive above 14 is preserved at [milestones/14/submission.zip](milestones/14/submission.zip).
