# Experiment 3: results, failed approaches, and the sign-mixing test

This report covers the three optimization tasks in experiment 3, including the
follow-up test of Claude's sign-mixing penalty. Some circuit builders were
inherited from experiment 2. Improvement factors below use the competition
baselines, rather than claiming that every improvement originated in this directory.

## Current result

The ready-to-upload [submission.zip](submission.zip) has a fully verified overall
score of **9.00962290835304**, up from **9.008540688342578** before this follow-up.
The inverse-square-root and DQI files are unchanged. The chemistry factorization
and its one-body angle-loading circuits were rebuilt and verified.

| Task | Competition baseline cost | Current verified cost | Improvement | Full verifier shots |
| --- | ---: | ---: | ---: | ---: |
| Inverse square root | 2,076,774.233715 | 190,711.955455 | 10.889586× | 200,540,928 |
| DQI | 24,439,463,583,672.93 | 1,612,737,419,463.0198 | 15.154025× | 102,400 |
| DFTHC-SOS chemistry | 3,355,599,588.315537 | 757,163,489.498628 | 4.431803× | 50,880 |

Lower task costs are better. The overall score is the geometric mean of the
three baseline-to-submission improvement factors. The task objectives differ:

- Inverse square root: scored qubits times Toffoli-equivalent cost, including
  rotation synthesis at the remaining approximation-error margin.
- DQI: repetition factor times qubits squared times Toffolis.
- Chemistry: effective normalization, λ_eff, times Toffolis times scored qubits.

Consequently, fewer gates alone did not reliably predict the best submission.
The exact sources and archive checksum are in
[submission.summary.json](submission.summary.json).

## Follow-up: does the mixing penalty help?

**It helped slightly at a weaker strength, but did not produce the predicted
λ_eff ≈ 48.4 or overall score ≈ 9.04.** The selected candidate has
λ_eff = **48.895187318**. The matched unpenalized refit reached **48.895751732**;
the previous submission had **48.912811162**.

The chemistry cost improved **0.036031%** over the previous submission. Only
**0.001154%** is the additional improvement over the matched unpenalized
control. Most of the gain came from continuing the refit. The overall increase
is **0.001082220**; the increase attributable to beating that matched control
is approximately **0.000034667**.

### What was measured

I used the penalty from
[Claude's fit_mix32.py](../3b/chem/fit_mix32.py), read without modifying experiment 3b:

`M(w) = Σ_c [‖w_c‖₁² − (Σ_b w_bc)²]`.

For a column with total positive mass p and total negative magnitude n, its
contribution is 4pn. The penalty therefore measures cancellation between signs.
It is invariant under flipping an entire column's sign.

The starting measurements reproduced the proposed comparison:

| Factorization | Shape B/C | Angle precision | Mixing M | Σ‖w_c‖₁² |
| --- | ---: | ---: | ---: | ---: |
| Previous experiment-3 submission | 719/350 | 8 bits | 81.08235 | 603.02953 |
| Claude's lowmix511_off array | 511/375 | 11 bits | 41.67816 | 563.63622 |

The second row changes the dictionary, column count, and angle precision. It is
useful motivation, but does not isolate what a penalty can achieve on the fixed
eight-bit dictionary. The measurements are saved in
[mixing_start_measurements.json](results/mixing_start_measurements.json).

### Controlled experiments

All refits started from the same verified `residual719_lbfgs2_onebody.npz`.
They held U, the Givens tree, all two-body angle words, B, C, and bit precision
fixed. Only the weights and identity offsets were optimized. Every run used
16,000 Adam steps, initial learning rate 2×10⁻⁶, whitened factor coordinates,
the w_B offset parameterization, and a two-body error target of 0.019995 Ha.
Checkpoints had to meet the physical error and SOS-gap acceptance margins.

Coefficient scaling matters. Claude's original objective is a physical squared
tensor error plus μM. The existing constrained refit uses λ_eff²/3000 and
normalized augmented-Lagrangian constraints. I tested both:

- **Error-scaled coefficients:** add `(μ / cap²) M` to that normalized objective,
  for μ = 10⁻⁷, 10⁻⁶, and 3×10⁻⁶.
- **Literal coefficient:** add `10⁻⁶ M` directly to the existing objective.

These are explicitly different conventions; μ is not assumed to transfer
numerically between the two programs. A zero-penalty run supplied the matched
control. The coefficient and effective scaling are recorded in each log.

I also ran matched 8,000-step orthogonal-column-rotation tests, with and without
mixing regularization. Rotating both w and its offset preserves the represented
two-body tensor, one-body correction, and SOS gap. This tests the penalty without
spending any physical-error budget. The regularized rotation coefficient was
0.9578 relative to the H2 normalization objective, approximately matching the
error-scaled μ=10⁻⁶ refit at the starting gap.

| Experiment | Checkpoint shown | Mixing | Verified λ_eff | Verified chemistry cost |
| --- | --- | ---: | ---: | ---: |
| [Previous submission](results/residual719_lbfgs2_verify.json) | starting fit | 81.08235 | 48.912811162 | 757,436,402.488 |
| [Refit: no mixing penalty](results/mixing_control_verify.json) | best cost | 80.20480 | 48.895751732 | 757,172,229.694 |
| [Refit: μ=10⁻⁷, error-scaled](results/mixing_mu1e7_verify.json) | best cost | 79.77083 | 48.895187318 | 757,163,489.499 |
| [Refit: literal coefficient 10⁻⁶](results/mixing_literal_mu1e6_verify.json) | best cost | 80.20975 | 48.895664221 | 757,170,874.554 |
| [Refit: μ=10⁻⁶, error-scaled](results/mixing_mu1e6_lowmix_verify.json) | lowest feasible mixing | 77.26907 | 49.016265447 | 759,038,437.602 |
| [Refit: μ=3×10⁻⁶, error-scaled](results/mixing_mu3e6_lowmix_verify.json) | lowest feasible mixing | 73.60242 | 49.418679944 | 765,270,003.147 |
| [Orthogonal rotation: no penalty](results/mixing_rotation_control_verify.json) | best H2 norm | 80.61609 | 48.903378503 | 757,290,333.609 |
| [Orthogonal rotation: mixing penalty](results/mixing_rotation_mu1e6_verify.json) | best H2 norm | 80.10585 | 48.912104861 | 757,425,465.101 |

All displayed rows have full native `PASS` results. The two stronger refit rows
show their **lowest-mixing feasible checkpoints**, rather than their best-cost
checkpoints. Neither stronger run beat its starting cost, so a best-cost selector
would retain the original fit. The weak, literal, and control rows show their
best-cost checkpoints. This distinction prevents a lower mixing statistic from
being mistaken for a submission improvement.

The selected weak-penalty candidate has mixing **79.77083**, one-body error
**11.996804 mHa**, two-body error **19.995290 mHa**, and SOS gap
**3.915138240 Ha**. Its Toffoli count and scored width remain **29,440** and
**526**, so its cost improvement comes from normalization alone.

### Why the larger reduction did not materialize

Mixing is a useful surrogate, but the score also depends on the offset terms,
the one-body normalization, and the SOS gap:

`λ_eff = sqrt(E_gap × (2Λ − E_gap))`.

At error-scaled μ=10⁻⁶, mixing fell to **77.26907**, but the gap rose from
approximately **3.91515** to **3.93929 Ha**. The verified λ_eff increased to
**49.016265447**. At 3×10⁻⁶, mixing fell further to **73.60242**, while the gap
rose to **4.01665 Ha** and λ_eff increased to **49.418679944**. The gain in raw
normalization did not compensate for the gap increase.

The rotation test reached the same practical conclusion without changing the
gap: its regularized candidate lowered mixing more than its control, but lowered
the full H2 normalization less. The literal-coefficient test also failed to beat
the selected weak error-scaled penalty.

The new submission uses the weakest error-scaled penalty. These are local results
for one warm start and fixed experiment budgets, not a proof of global
optimality or a newly established lower bound on λ_eff.

Machine-readable parameters, all measurements, verifier paths, and comparison
scores are in [mixing_experiment_summary.json](results/mixing_experiment_summary.json).
The implementations are [constrained_factor.py](constrained_factor.py),
[rotate_qr.py](rotate_qr.py), and [measure_mixing.py](measure_mixing.py).

## Task 1: inverse square root

### Ideas that worked

**Binary normalization and centered polynomial interpolation.** Normalizing the
mantissa and folding normalization parity into the coefficient lookup reduced
the interval count. Searching nearby integer coefficients, using only reachable
normalized inputs, and retaining synthesis-error margin produced cheaper valid
approximations than simply increasing coefficient precision.

**Short accumulators with delayed carry integration.** Horner accumulation used
the low portion of each word plus a small overflow counter. Its signed carry
was integrated into the high portion once per stage. Rounded partial products
allowed narrower coefficients while preserving the required phase accuracy.

**Radix-four Booth accumulation.** Sharing a ripple addition between adjacent
coordinate bits lowered arithmetic cost. The implementation restores its controls
and erases temporary products with the required phase corrections.

**Trading gates for workspace.** Hybrid adders use in-place carry propagation
where necessary to fit a qubit budget, and temporary-AND carries elsewhere.
Borrowing idle output bits, reusing the normalized input's known top bit, and
delaying allocations improved the scored qubit–gate product.

**Delayed coefficient phase repair.** After the first Horner stage, unused high
curvature bits are measured and released. Their measurement phase is repaired
when the coefficient is reconstructed for the inverse stage, using conditional
Z gates rather than a separate phase lookup. Carry borrowing and private-wire
lifetime compaction then reduce the circuit to **94 qubits**.

| Full verified stage | Cost |
| --- | ---: |
| [Binary-centered interpolation](results/binarycentered_verify.json) | 298,066.978 |
| [Rounded partial products](results/rounded_f8_w27_16_8_verify.json) | 233,305.885 |
| [Adaptive Booth table](results/adaptivebooth_f7_m8_w27_18_8_verify.json) | 209,396.519 |
| [98-qubit workspace schedule](results/budgetunusedbooth_q98_verify.json) | 194,025.301 |
| [Delayed phase repair and compaction](results/compactdefer_k4_q94_verify.json) | 190,711.955 |

The final circuit uses approximately **1,778 mean Toffolis**, **2,028.851
Toffoli equivalents** after synthesis, and 34 rotations. Exhaustive small-adder
tests checked values, phases, and workspace cleanup; the final artifact passed
all **200,540,928** native shots. [build_compact_invsqrt.py](build_compact_invsqrt.py)
reproduces its bytes.

### Ideas that did not help

- **Dropping a low phase rotation:** passed verification, but consumed
  approximation margin and increased the remaining rotations' synthesis cost.
  [Result](results/invsqrt_drop0.json).
- **Global phase offsets, alternate tie rules, and rounding-mask searches:**
  did not yield a cheaper passing candidate in the tested grids. Reducing
  coefficient precision often left a small set of intervals outside tolerance.
- **Immediately repairing an erased curvature coefficient:** the extra lookup
  cost offset the saved qubits. Delaying repair was the useful change.
- **Boolean synthesis of the high curvature bits:** the tested shared expansion
  still required 41 ANDs and was not competitive with the chosen lookup schedule.
- **Minimizing Toffolis alone:** some wider circuits had fewer arithmetic gates
  but a worse final score. Live workspace and synthesis margin both mattered.

The coefficient grids and failed variants remain in `results/`; the sequence of
changes is documented in [CONTINUATION.md](CONTINUATION.md).

## Task 2: DQI

### Ideas that worked

**Adapting a public decoder rather than rebuilding every primitive.** The
attributed upstream implementation in [vendor/uncompute](vendor/uncompute)
provided a reversible Berlekamp–Massey framework. The adaptation supplies runtime
masks and reversible conversion between the public and host field bases. The
upstream license and revision are retained.

**Consume syndromes, then reuse their storage.** Once a syndrome is no longer
needed in its original form, the decoder consumes it and uses its register for
the backup polynomial bank. A parked input handles a vanished leading
coefficient; reversible stop guards protect exceptional paths. Scheduling bank
materialization after consumption reduced the peak allocation.

**Recursive finite-field phase operations.** Subfield norm and inverse-parity
contractions made the Forney phase scan cheaper. The recursive kernel costs
80 Toffolis per evaluation point. An implicit monic coefficient avoids storing
a separate constant word.

**Select the six-bit inverse norm before the final product.** This avoids
materializing separate twelve-bit inverse and intermediate words. Using the
armed locator coefficient directly as a scale removes another duplicate register.
Together these changes brought the final width to **1,429 qubits**.

| Full default-size verified candidate | Cost |
| --- | ---: |
| [Initial published adaptation](results/published_verify.json) | 1,829,622,334,585.752 |
| [Syndrome-consuming stream3](results/stream3_full_verify.json) | 1,772,681,819,628.758 |
| [Corrected final-round guards](results/stream4_full_verify.json) | 1,681,616,324,033.482 |
| [Final compact selection/scale schedule](results/stream10_offset3_select_noscale_compact_verify.json) | 1,612,737,419,463.020 |

The final circuit uses **645,755 Toffolis** and a scored repetition factor of
**1.223013998**. Its 102,400-shot result has **zero value errors and eight phase
errors**. The host accepts and prices those phase errors through its repetition
model; `PASS` does not mean that this decoder had zero empirical phase errors.

### Ideas that failed or needed correction

- **Early streaming implementations:** journal replay, iterator lifetime, and
  translation bugs broke uncomputation. In one case an empty-file insertion
  stopped output and omitted the reverse pass. These were implementation
  failures, fixed before selecting a submission.
- **A stale final-round guard:** incorrectly stopped valid weight-35 cases.
  Tracing failing shots identified the obsolete sign/stack check; removing it
  reduced the stream3 failures substantially.
- **A smaller bank without extending syndrome consumption:** produced roughly
  1,400 phase errors. Gating discrepancies added about 32,000 Toffolis and did
  not resolve the cause. Matching consumption to the actual bank lifetime did.
- **Compaction without compatible replay scheduling:** could create allocation
  conflicts. The final construction prematerializes the needed locator storage
  before same-round scratch is freed.
- **Extra cofactor workspace and alternate field formulas:** no tested option
  beat the final total cost; extra registers could erase an arithmetic saving
  because width is squared in this score.
- **Small-instance apparent wins:** tests using fewer evaluation points were
  useful diagnostics, but their scores were not substituted for the full
  4,095-point task.

Correctness diagnostics and resource accounting were both necessary: some
experiments returned native `PASS` yet remained poor submissions because their
error rate or workspace increased the priced cost.

## Task 3: DFTHC-SOS chemistry

### Ideas that worked

**Compressing the inherited ensemble factorization.** Reducing column counts
within the existing 35-ensemble representation produced early valid savings.
The C=50 candidate cost about 1.773 billion, and further constrained refinement
reached about 1.758 billion.

**One shared dictionary with PSD projection in the correct metric.** A global
dictionary substantially reduced the number of angle rows. The physical tensor
loss lives in a symmetric feature space; projecting the coefficient matrix in
ordinary Euclidean coordinates was inappropriate. PSD fitting in the feature
Gram metric, followed by evaluation against the actual reference tensor, made
this route usable.

**Optimize the lattice, not just the floating-point fit.** Integer-angle moves,
rank-constrained refits, and growing or reseeding the quantized dictionary
eventually produced a valid eight-bit dictionary with 719 vectors. Floating-point
accuracy alone was an unreliable guide to quantized accuracy.

**Optimize normalization while preserving or constraining physics.** Physical
column alignment and orthogonal rotations reduce LCU cost without changing the
represented tensor. Identity-offset optimization lowers the SOS gap; a small
PSD update can spend remaining tensor-error budget to approach its acceptance
floor. Subsequent constrained weight/offset refits, especially in whitened
coordinates, reduced the effective normalization further.

**Optimize one-body angle words separately.** Direct rounding was insufficient
for some otherwise useful candidates. Joint integer-row optimization brought the
one-body error below its threshold. After changing weights or offsets, these
words were regenerated; exact column rotations could reuse them because their
one-body matrix is invariant.

**Price the complete preparation circuit.** Eight-bit precision changes the
state-preparation resource regime. Six batches of thirteen angle words fit the
128-qubit per-circuit cap. The final preparation cost is 29,440 Toffolis and
526 scored qubits, even though an individual angle-loading circuit is narrower.

| Full verified stage | Chemistry cost |
| --- | ---: |
| [35 ensembles, C=50](results/chem35c50_verify.json) | 1,773,248,191.824 |
| [Global dictionary, 11-bit polished fit](results/global500_factor_polish_verify.json) | 909,956,199.771 |
| [10-bit fit](results/lbfgs500_q10c375_verify.json) | 838,859,431.309 |
| [Initial valid eight-bit 719-vector fit](results/residual719_q8_initial_verify.json) | 823,241,449.185 |
| [Previous constrained refinement](results/residual719_lbfgs2_verify.json) | 757,436,402.488 |
| [Weak mixing-penalty follow-up](results/mixing_mu1e7_verify.json) | 757,163,489.499 |

### Ideas that did not produce a better submission

- **Aggressively reducing ensembles, basis size, or rank:** some fits passed in
  floating point but failed after quantization. Other low-rank fits stalled above
  20 mHa. In a late 350-column snapshot, removing one physical column gave
  approximately 20.0095 mHa; keeping only 340 gave 20.330 mHa.
- **Adding columns or overcomplete factors:** small normalization reductions
  generally did not repay the extra state-preparation cost.
- **Straight-through quantization and unrestricted angle updates:** several
  runs were unstable or improved only their soft angles. A joint continuation
  ended near 20.545 mHa with hard angles, still outside the gate.
- **Strictly nonnegative weights:** attractive sign statistics did not suffice.
  A positive 719-vector bound-constrained fit ended near 22.064 mHa; an enlarged
  988-vector version reached about 21.662 mHa. Both remained invalid.
- **Multi-ensemble block fits:** a floating fit could reach about 14.068 mHa,
  yet its eight-bit hard fit remained around 31 mHa. Splitting the final global
  factors into smaller clustered blocks produced much larger projection errors.
- **Direct use of published THC factors:** their raw spans did not fit this
  benchmark's fixed symmetry shift accurately enough. They were useful as new
  initializations, rather than immediate submissions. Attribution is retained
  in [vendor/published_thc/SOURCE.md](vendor/published_thc/SOURCE.md).
- **Large random rotations, proximal rotation schemes, tree proxies, and tiny
  atom swaps:** none produced a better verified final candidate in the tested
  configurations. A better geometric proxy did not guarantee a better physical
  tensor fit.
- **Stronger mixing penalties:** this follow-up confirms that lower cancellation
  alone can worsen the actual objective. The weak penalty gave only a small
  local benefit over a matched continuation.

Some early apparent improvements were also numerical or implementation errors:
naive PSD clipping changed a roughly 19 mHa fit to roughly 68 mHa, and an
incorrect spectral ordering invalidated one accelerated fitting experiment.
Those artifacts were excluded. Final acceptance always used the actual reference
and the unchanged host verifier.

## Verification, provenance, and limits

Seven new chemistry candidates in the mixing study each passed the full
50,880-shot verifier. The unchanged inverse and DQI components reuse their
existing full-size `PASS` results. Array hashes confirm that all mixing-study
candidates retain the same U, tree, and two-body angle words. Only weights,
offsets, and where necessary regenerated one-body words changed.

The [final archive audit](results/final_archive_audit.json) binds all nine archive
members to their verified source files and checks the chemistry submission digest
and preserved milestone checksums. Milestones **5–9** remain intact; milestone 9
is the original 9.008540688 checkpoint, while the current archive contains this
later polish. A copy of the new archive is retained under
[revisions/mixing_mu1e7](revisions/mixing_mu1e7).

The evidence supports the reported artifacts and tested local improvements. It
does not establish global optimality, prove every unsuccessful family impossible,
or turn short diagnostic runs into full-task verification. Complete historical
logs remain in `results/` and [CONTINUATION.md](CONTINUATION.md).
