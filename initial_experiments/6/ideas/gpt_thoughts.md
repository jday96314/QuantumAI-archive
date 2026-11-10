# Quantum Circuit Compilation: Autonomous Research Backlog and Breakthrough Playbook

**Prepared:** October 24, 2026  
**Purpose:** Read-only research agenda for multiple long-running `/goal` agents attempting local combined scores of **13, 17, or 20+** in *Google Quantum AI — Verifiable Quantum Circuit Compilation*.  
**Scope:** All three tasks (`dqi`, `dfthcsos`, `invsqrt`), their integration, and the search process itself.  
**Interpretation:** A large, deliberately diverse collection of **hypotheses to test**, not a claim that any of the proposed improvements is achievable. Source reports are listed in Appendix A. Public papers are linked in Appendix B.

> **Central instruction to autonomous agents:** Do not convert a failed implementation, an optimizer plateau, or a restrictive choice of internal coordinates into a claim that the competition target is impossible. Equally, do not report invented improvements or treat an unbuilt projection as a verified circuit. When an idea fails, determine precisely *which representation, operator, accuracy constraint, or resource term failed*, preserve the evidence, and search a materially different realization.

## 0. Immediate orientation: what is already known

### 0.1 Best *separately verified* task costs found in the supplied reports

| Task | Best reported cost (lower is better) | Source | Important distinction |
|---|---:|---|---|
| Inverse square root | **182,956.5271606425** | `initial_experiments/4/codex`, `invsqrt` selected circuit | 90 quantum wires; full 200,540,928 shots |
| DQI | **1,435,930,717,449.9397** | `initial_experiments/5/dqi-gap/candidates/gap_best_n70_m4095.kmb` | 1,004 quantum wires; no errors in three separately reported 102,400-shot full checks |
| Chemistry | **696,942,846.9115** | `initial_experiments/5/chem-census/candidates/structured` | 454 scored qubits, 31,608 total Toffolis, `lambda_eff=48.56733747`; full native chemistry verification twice |

Those costs imply a **calculated** local score near **9.76167123** when combined, using the same host denominators as the previous 9.36772215 checkpoint. **This combination was not shown in the supplied reports to have been assembled and checked in one complete metric run**. Verify integration before promoting it. A reported fully packaged archive scored **9.47124939**, but used the inferior (for the recorded default samples) `dqi-core/trimmed_decoder.kmb`, not the gap-aware decoder. Do not silently replace the best DQI component with the latest *packaged* DQI component.

The improved chemistry result comes from combining (i) the `chem-audit` physical normalization checkpoint (`lambda_eff=48.56733747`) with (ii) the `chem-qrom` four-row structured compiler / smaller batch size. The uncombined `chem-qrom` physical factorization scored 697,327,965.84. The earlier 708,156,073.79 and 708,547,388.96 chemistry scores remain useful controls.

### 0.2 Score arithmetic; why 13, 17 and 20 need qualitatively different strategies

Let `g_D`, `g_C`, and `g_I` be **additional** cost-reduction factors, each relative to the best costs above. Then

```
projected_local_score = 9.76167123 * (g_D * g_C * g_I) ** (1/3)
```

| Goal | Required additional product `g_D*g_C*g_I` | Illustrative factor combination (DQI × chemistry × inverse) | Projected score |
|---:|---:|---|---:|
| **13** | **2.362×** | 1.6 × 1.35 × 1.1 | ~13.03 |
| **17** | **5.282×** | 2.5 × 1.8 × 1.2 | ~17.13 |
| **18.9** (user-reported leader) | **7.258×** | 3.0 × 2.0 × 1.2 | ~18.85 |
| **20** | **8.600×** | 3.5 × 2.1 × 1.2 | ~20.17 |
| **20** | **8.600×** | 3.0 × 2.0 × 1.5 | ~20.31 |

These are **arithmetic illustrations, not predictions**. The public leader's reported ~18.9 is evidence that better complete scores exist, but does not reveal whether the gains are in DQI, chemistry, inverse, or some combination. Do not assume every task must contribute equally. In particular, do not infer an upper bound on the combined score from the exhaustion of one local architecture.

**Suggested allocation:** For a target of 13, favor opportunities with a concrete path to complete native verification; for 17, require at least one architecture-scale gain; for 20, prioritize distinctly different algorithms and new representations alongside banking cheaper wins. Agents may pursue a different path if their evidence supports it.

### 0.3 Host contracts and genuine constraints

- `dqi`: RS(4095,4025), 70 quantum 12-bit syndromes (**840 unavoidable input qubits**), 4,095 **classical** 12-bit measurement masks, arbitrary per-shot masks; correct syndrome restoration, clean scratch, appropriate parity phase. Score `R(epsilon_plus)*Q**2*mean_Toffolis`; host uses 102,400 shots and one-sided Clopper–Pearson with `delta=2^-128`. Even zero observed failures gives `R≈1.19126046`. The public Gu circuit's 442k Toffoli / 1,837-qubit figure cannot be transplanted naively because this host squares the qubit count.
- `dfthcsos`: 76 spatial orbitals; physical vectors **decoded from Givens words** rather than assumed equal to an NPZ's floating-point `U`. The NPZ specifies `U(R,B,76)`, `W(R,B,C)`, `H2_I(R,C)`, one-body fields, and a shared 75-edge tree. Both physical residuals must be <=0.020 Ha, with strict `Lambda > 49`, `E_gap > 3.914898` and positive allowed domain; fixed native symmetry shifts. ANGLE family must implement cumulative XOR *differences*, arbitrary incoming targets, clean scratch, and each oracle <=128 qubits. Score `lambda_eff*T_total*Q_total`; `T_total` includes actual ANGLE and modeled PREPARE/Givens; `Q_total` is persistent plus the **maximum shared temporary** footprint, not a sum of every oracle's width.
- `invsqrt`: every integer `x=1..3,133,452`, 64 native measurement realizations each (**200,540,928 shots**), target phase proportional to `1/sqrt(x)`, error tolerance `2^-24/sqrt(x)` turns, input restored and scratch clean. Score `Q*(mean_Toffolis + total_T/4)` with synthesis costs depending on the **smallest absolute remaining phase margin**, not just the worst relative error. One extra rotation can cost more than it appears after remargining.
- The **universal** native chemistry normalization-only bound, for this fixed shift and SOS model, is `lambda_eff > 46.7861795925` (`chem_audit_assumptions_and_margin_report.md`). It leaves only ~3.86% normalization headroom from `lambda_eff=48.5673`; it says nothing about reducing `T*Q` or alternative legal factorizations. The fresh symmetric-pair spectrum also gives the numerical necessary **rank `R*C >= 290`** for a two-body residual <=20 mHa in this PSD factor model; passing that condition is far from sufficient.
- All calculations should use the **unchanged host implementation**. The canonical input/phase contracts and score matter; simulator runtime, conditional Clifford count and archived agent estimates are *not* the score, though runtime matters for exploration speed.

### 0.4 Lessons from the latest 15 reports

| Evidence | Conclusion to retain | Scope of that conclusion |
|---|---|---|
| `dqi-gap` | A coherent gap-aware decoder **works**: 1.436T cost, zero errors in three full native runs | Arithmetic now 1.196M mean T / 1,004Q: correctness alone was not a large enough gain |
| `dqi-core` / `DESIGN_MEMO.md` | Trimming zero middle-product rows and retaining 44 inverse chains saves exactly 25,596T at 961Q; trimmed candidate fully verified | Still has ~380 phase errors; not superior to gap-aware DQI at zero-error R |
| `dqi-dispatch` (both variants) | Fixed-schedule exact classical transfers work across gaps; some shared quantum components work | Separate moment charts / Bennett bridges add 176k–350k+ T, and width can exceed target. **Do not rerun those same implementations.** |
| `dqi-core-phase` | Validated Riccati `DS=S²`, `DB=Omega`, `D Omega=0` in the supported Padé frame; in-place C-fraction small decoders and a 46T inverse-only point body | Full C-fraction extraction costs ~1.123M T before scan, generic chart failures too common; 46T/H phase needs better upstream recurrence |
| `dqi-census` | 960T / +69Q **degree-two** inverse norm primitive; rational-coordinate and direct F-prime identities verified | Degree-35/odd-degree generalization missing; optimistic 2.65× DQI scenario is **a design target**, not a measured full pipeline |
| `chem-qrom` | Native-verified 1.016× chemistry improvement by moving J=7/Q478 to J=4/Q454 with an exact nonlinear compressed bank | Dense source table defeats simple repeated-subtree/bitwise-affine compression; higher-order bank decomposition remains open |
| `chem-census` | Compatible `chem-audit`+`chem-qrom` candidate native-verified at **696.943M** | Physically smaller B703 becomes valid only with much worse normalization; simple prune/grow does not preserve low-Lambda basin |
| `chem-audit` | Native-verified 708.156M normalization improvement; universal lower bound scoped correctly | Global normalization-only breakthroughs are mathematically blocked *under this verifier*; resource changes are not |
| `chem-tree` | Legal tree variants tested; modest quantization gains; near-singular directions destroyed by eight-bit rounding | Specific tree/word tables ruled out by span/PSD certificates; **new robust discrete geometries** still open |
| `chem-psd` / `chem-discrete` | Continuous R2/B351/C189 achieves ~19.98 mHa; fixed rounded dictionary's union-span residual >140 mHa; B479 has fixed-word PSD dual lower bounds >20 mHa | Do **not** mistake fixed-dictionary capacity bounds for shape-wide impossibility. C/rank relaxation cannot repair a bad span |
| `inv-dyadic` | Independently discovers almost exactly incumbent's 269-leaf prefix table; full native alternatives tie/lose | Another trivial interval split is unlikely to help; no lower bound on other arithmetic algorithms |
| `inv-phase` | Direct row phasing and split-bank phase circuits work but pay many rotations/shift gates | Numerically passing **exact rational quadratic** has ~14.8× incumbent synthesis margin; compact circuit for it is unresolved |
| `inv-census` | **Unbuilt** narrow fused-square correction passes exhaustive numeric test; injective exponent marker also proved numerically | This is one of the strongest unfinished inverse arithmetic leads; count full arithmetic and cleanup before predicting cost |
| `inv-width` | 90Q peak is inside the first row of last Horner product; 88Q builds lose on extra gates | In-place normalized input and in-place coefficient Horner already implemented; borrowed dirty-ancilla and ordinary split-add do not win in tested schedules |

---

## 1. Prioritized campaign map — not a universal ranking

**A: near-term / better measured success probability.** Pursue at least one of these in each long-running research wave:

1. **A1 — Assemble and validate best independent task components.** Use `dqi-gap`, `chem-census/structured`, and round-4 inverse; score the complete archive with unchanged full metric. Bank any valid improvement; keep a frozen fallback.
2. **A2 — Exploit the new chemistry J4 quadratic-bank compiler.** Jointly optimize which columns are omitted, their shared nonlinear reconstruction, low-address grouping, output-affine basis, row permutation, and batch size. The remaining compressed lookup columns are a *different* search space from repeated-subtree pruning.
3. **A3 — Inverse fused-square correction.** The numerical result in `inv-census` is not a built circuit. Prototype the exact narrow square, final multiply/phase, and live workspace; compare to the incumbent rather than to a generic squarer.
4. **A4 — Reuse DQI middle-product trimming and phase-mask interpolation wherever compatible.** Do **not** assert trim and gap code compose unchanged: reproduce the exact code and verify both branching and resource lifetimes.
5. **A5 — Multi-stage cost frontier and component-level auto-tuning.** For each task, find physically admissible or classically correctness-tested implementations, and use exact host resources to choose winners rather than trusting hand estimates.

**B: plausible *architectural* breakthroughs.** These should dominate the 17/20 research budget:

6. **B1 — In-register DQI locator/evaluator or C-fraction representation**, avoiding duplicated polynomial banks *and* repeated whole-series reciprocals.
7. **B2 — Joint DQI phase-aware recurrence** that evolves an object equivalent to F, F-prime or the Wronskian H *during decoding*, not via a late wide converter.
8. **B3 — Cost-optimized approximate DQI decoding**, deliberately trading rare degree-gap errors for simpler arithmetic; sweep error **incidence by stage** and host `R`, not arbitrary thresholds.
9. **B4 — Robust quantized chemistry from inception:** new dictionaries with well-conditioned high-leverage spans and low LCU normalization, rather than rounding fragile rank-deficient float solutions.
10. **B5 — Rank/C/B/R shapes that cross host width cliffs while preserving chemistry accuracy**, including low-rank shared dictionaries and compressed physical oracles.
11. **B6 — Phase-aware inverse-square-root arithmetic:** fused exact rational quadratic, narrow partial products, shared exponent routing, and stage-wise phase application with a **small number** of rotations.

**C: high-risk, high-upside explorations.** Give these time for theoretical development and small-case experiments, not immediate full-circuit builds:

12. **C1 — Reversible Padé/half-GCD/partial-inverse with branch history reconstructed from state.** Seek one inexpensive algebraic invariant; do not port a wide classical implementation wholesale.
13. **C2 — Phase as a syndrome/mask contraction**, using compressed polynomial remainders or transposed operations; avoid scanning all 4,095 points if total quantum work permits.
14. **C3 — Joint chemistry factorization/lookup-code co-design**, paying a small physical or normalization penalty for a dramatic reduction in Boolean lookup complexity.
15. **C4 — Alternative inverse basis or digit-recoding** that makes the phase a cheap sum of few *charged* rotations rather than 201 independently shifted terms.
16. **C5 — New reversible finite-field bilinear circuits and specialized matrix/series kernels** tailored to the actual operand shapes and scratch windows, not generic scalar multipliers.

For each item below, the labels mean: **E** = grounded by measured evidence; **P** = partially built; **S** = speculative. They rank practical usefulness under current evidence, not mathematical possibility. Every idea has a recommended *first test* and a clear warning about what has already failed.

## 2. DQI: 27 research directions

**Current anchor:** gap-aware `Q=1004`, mean `T≈1,195,802`, zero observed failures, `R=1.19126046`, cost **1.43593e12**. The arithmetic core is ~824,740 static Toffolis, scan ~371,150. The older 961Q trimmed decoder has ~994,270 executed Toffolis but ~380 phase errors and is therefore worse under the confidence penalty in the reported full runs. Use *both* as Pareto controls.

**DQI break-even worksheet at `R=1.19126046` (zero failures at default shots):**

| Q | Mean Toffolis for 1.3× vs gap decoder | For 2× | For 3× |
|---:|---:|---:|---:|
| 840 | 1,314,089 | 854,158 | 569,439 |
| 900 | 1,144,718 | 744,066 | 496,044 |
| 940 | 1,049,368 | 682,089 | 454,726 |
| 961 | 1,004,007 | 652,604 | 435,070 |
| 1,004 | 919,848 | 597,901 | 398,601 |
| 1,050 | 841,017 | 546,661 | 364,441 |

At nonzero phase errors, recompute `R` from the host Clopper–Pearson function. Do not use zero-error `R` for projected architectures with uncorrected generic exceptions.

### D01 — Blocked in-register C-fraction coefficient extraction [B1, P, potentially transformative]

**Rationale:** `dqi-core-phase` built an in-place C-fraction decoder that uses the original 840-bit syndrome bank plus modest overhead and produces correct small-instance locator/evaluator pairs on supported charts. Its serial N70 extraction alone costs ~503,655 T forward, with ~57,664 T for synthesis; extraction+reversal is ~1.123M T before the phase scan. The *space representation* is promising; the per-coefficient reciprocal algorithm is not.

**First test:** Derive a 2/4/8-coefficient **blocked C-fraction/Padé step** using one short reciprocal or a shared matrix product, whose outputs live in exactly the consumed syndrome words. Implement N6/N14 and N22 controls first; benchmark gates per emitted coefficient and transient live words. Require correctness through zero-head/degree-gap charts on host-law samples, or explicitly price residual failures. Compare against the serial baseline at equal width and reversible cleanup.

**Avoid:** Repeated whole-70-series reciprocal per coefficient, separate sigma/template/tau banks, and the assumption that a small-instance `Q=878` profile automatically scales with no extra registers. The gain requires a **different blocked arithmetic recurrence**, not just another circuit generator for the same recurrence.

### D02 — Phase-aware locator/evaluator Riccati tangent [B2, P, mathematically intriguing]

`dqi-core-phase` independently verified for supported exact Padé pairs that the syndrome series' derivative identity `D S = S²` leads to tangent coordinates `D B = Omega`, `D Omega = 0` (check conventions and characteristic two carefully). This suggests an evaluator-containing *differential* decoder state rather than synthesizing Omega or `F` after decoding.

**First test:** Symbolically derive the derivative of **one actual quotient/Schur update** and determine whether the tangent can be stored in an already-resident cofactor/head coordinate. Verify exact finite-field algebra on 0–35 errors with degree-two and degree-three quotient cases. Emit a small reversible tangent-update circuit; tabulate Toffoli cost per stage, one-pass restoration and additional **peak** words. If keeping a full derivative polynomial is required, find a transposed/adjoint evaluation that transports only the phase-relevant contraction.

### D03 — Jointly update locator B and Wronskian H [B2, P]

`dqi-core-phase` produces a degree-35 H preparation in **100,080T build/erase at ~897Q including inputs**, plus a verified **46T** inverse-only scan body. A naive completed adaptation has a very small zero-error margin and generic-chart failures. Rather than preparing H after finishing B, maintain `(B,H)` directly through quotient updates, using Wronskian/Sylvester identities and reconstructing any unneeded A information on the way back.

**First test:** Derive the exact `H -> H'` update under one generic degree-one and degree-two Euclidean step without an expensive inverse Sylvester reconstruction of A. Prove reversible recovery of the overwritten coordinate, then measure a complete 2–3 stage component including cleanup. **Avoid** the prior separate H-bank scheme and repeating the 50,040T/direction postprocessing unchanged.

### D04 — Direct F-prime evolution with classical mask factors [B2, P]

At locator roots, `e(beta)=beta^70/[A(beta)B'(beta)]`. Compute `F'=(A B')^-1 mod B` and absorb `beta^70` into the *classical* mask's trace-dual factor. `dqi-census` verified this identity across all degrees and masks. Its degree-two inverse norm costs **960T and +69Q**, but straightforward Frobenius inversion at degree 35 is disastrously wide (11,751Q).

**First test:** Search for a split-algebra/Bezout/Hensel **in-place** update of `F'` with sublinear scratch. A valid degree-2 norm formula is a base case, **not** a degree-35 algorithm. Prefer a phase-aware recurrence that avoids computing `F'` by a late quantum-modulus inverse; profile the full quotient-history and scratch. Reject only the concrete broad Frobenius workbank, not modular inverses generally.

### D05 — Unified gap-aware moment arithmetic, not separate quantum branches [A/B, E+P]

`dqi-gap` achieved correctness but added arithmetic: degree-10 ordinary three-slot 34,953T versus gap-aware 43,470T; the new core total is ~824,740T. `dqi-dispatch`'s fixed classical decoder repairs degree-1/2/3 patterns, but its measured **separate moment charts** add at least ~326,622 T before all cleanup. A **single shared moment map** across ordinary, carried quadratic and cubic charts remains unbuilt.

**First test:** Factor all branch formulas into common symmetric-polynomial/Schur moments. Use conditional *linear* coordinate permutations to supply one uncontrolled nonlinear product path. Benchmark one normal/degree-2/degree-3 block at k=1,10,19,33, including both directions, and demand that routing overhead stay below the current gap-aware block. Do not implement the losing matrix/Bennett bridge again.

### D06 — Degree-gap tolerance as an explicit optimization variable [A/B, new strategic direction]

The gap decoder pays for essentially exact correction, even though DQI's host permits errors and prices them via `R`. The old compact trimmed decoder tolerates hundreds of phase errors for cheaper arithmetic. The best design might be **partially corrected**: maintain cheap ordinary updates and fix only the stages/cases with the highest *failure-reduction per marginal Toffoli/qubit*.

**First test:** Produce a stage × quotient-pattern census with incidence and mask-phase probability from >=2M independent law samples. For every correction subset, compile or conservatively count the quantum-selected gates (not occurrence-weighted arithmetic), predict `R` for 102,400 shots using the exact confidence bound, and optimize `R*Q²*T` by dynamic programming/knapsack plus Pareto search. Investigate correcting only early, late, or easy-to-share gaps; preserving rare degree-three failures; allowing a few dozen phase errors if it removes a large circuit block. **Caveat:** quantum controls still execute: savings require *deleting circuitry*, not conditioning it on a rare quantum event.

### D07 — Recode completion history from 35 unary flags into binary counter [A/B, P]

`dqi-dispatch` verified a reversible valid-completion-suffix counter: `22 flags -> 6 count + 3 local flags`, suggesting a ~13Q saving near the old block peak, and an estimated 959Q scan -> ~936Q with an explicit leading word. The original approximate decoder's spurious gap flags **violate the suffix promise**, so it is not a drop-in patch.

**First test:** Apply to the *correct* gap-aware decoder, for which the completion invariant differs. Prove its history tags and pending-carry modes support reversible count regeneration on all charts. Re-run full liveness and scalar phase scan, including the explicit monic word, and compare the savings after counter extraction/cleanup. Do not simply substitute six bits for flags while losing history.

### D08 — Merge the gap-aware decoder with middle-product trimming [A, E+P]

The 25,596T saving of `dqi-core` is fully verified on the ordinary compact decoder. The gap-aware decoder has different block layouts and extra carry forms. Re-derive coefficient-independent structural zeros for **each gap-aware block chart**; prune evaluated/transposed middle-product rows/columns *before* selecting a multiplication algorithm. Profile whether the `budget_inverse.py` chain fits in the 1004Q gap-aware liveness window.

**First test:** Component-wise exact replay and 256/512-shot checks at k=1,10,19,33, then full circuit. If a subset transfers, bank it even if the total saving is only 10k–30k T; don't assume the entire 25.6k transfers unchanged.

### D09 — Global degree-gap branch reconstruction, stack compression, and history-free inverses [B, P]

The public `uncompute` decoder reconstructs branch information from invariants instead of storing worst-case history. The gap-aware circuit still carries completion bits, carry data and stage-specific transformations. Identify variables that can be reconstructed on reverse traversal from the final cofactor pair plus one compact tag. Compare *stored decision*, *recomputed decision*, and *measurement-erased* variants on a per-stage basis. A proof of injectivity on all promised quantum-state fibers is required; zero-pivot overwrite can be many-to-one, as `dqi-hybrid` demonstrated.

### D10 — Save cofactor passes with Laurent/prefix arithmetic [B, E+P]

`dqi-gap/laurent_cofactor.py` removes conditional current-bank rotations. Its one-pass late-stage update reconstructs consumed remainder words from cofactor/normalizer invariants and avoids copied-head passes. Applying it early repeats long-tail work. Search for a **partial one-pass** that shares the repeated tail evaluation with adjacent blocks, caches only a few minimal coefficients, or uses a transposed product across two time steps. Measure complete forward and inverse arithmetic at early/mid/late k. Do not repeat the already-costed full ratio recurrence, which lost by ~3–10× early.

### D11 — Truly in-place (sigma, Omega), (sigma,F), or (sigma,F-prime) final coordinates [B, P]

`DESIGN_MEMO.md` gives an endpoint identity `P B + R Bprev = p x^70`. Reversibly trading `(Bprev,pinv)` for an evaluator uses **existing** 35-word storage in principle, but the converter's lifetime must not temporarily hold two 35-word banks. A modest 80T phase kernel saves only ~32,760T/scan and permits **<16,380T/direction** for its adapter in the old layout. The 11T F scan buys more budget but its known converters are too wide. Work backwards from an alias-safe, reversible state layout; prove erasure algebraically before producing a giant converter.

### D12 — General-quotient Euclid in a fixed register frame [B/C, E+P]

The two dispatch reports contain **classically correct** local degree-1/2/3 quotient transfer rules and clean component formulas. Their cost explosions came from translating six-moment coordinates to explicit 2×2 polynomial matrices and back. Keep the quotient *factored* (scalar/linear/quadratic coefficients and a small carry descriptor) and apply its action via one set of cofactor-tail MACs. Avoid expanding intermediate matrix entries; benchmark sparse factorized 2×2 actions against the measured 8,638T/18-pair matrix interior and 4,356T Bennett bridge.

### D13 — Sliding-block/half-GCD across 70 syndrome moments [B/C, S]

Degree 70 is small enough that asymptotically fast half-GCD may lose after reversal. But **subquadratic short middle products** or blocks of unequal size may win if they reuse right-hand polynomial evaluations and avoid a copied 6/8-word head. Build a costed local dynamic program over block widths 1..8 with *actual arithmetic/liveness* on generic and defective paths. Existing block4 peaked at 984Q against the old 961Q cap, but a different head layout can change that; don't reuse the same block4 schedule without a storage innovation.

### D14 — Reversible polynomial approximant-basis / Padé with minimal state [C, S]

Look at minimal approximant basis, partial-inverse, rational reconstruction, and Euclid/Berlekamp–Massey equivalences. Instead of forming a complete basis and later discarding one column, derive a **streaming 70-word representation** that carries only the needed transform plus enough invertible history. For each small-instance candidate, count expensive finite-field multiplications and required quantum words. A classical O(n log n) operation count is not automatically favorable at n=70; require quantum reversible resource accounting.

### D15 — Hankel/Toeplitz displacement-rank decoding [C, S]

Syndromes define structured Hankel matrices. A full 35×35 quantum linear-system solve is hopelessly wide; `dqi-census` found a 2200Q displacement realization and dense PGZ could need ~14,700Q. The open avenue is **generator-only rank revealing**: carry O(n) displacement generators rather than storing the matrix, derive cheap exceptional-pivot history, and extract a locator/evaluator from them. Reject only measured bulky realizations. Test compressed structured factorization at degree 7/15 before N70.

### D16 — Use transposed polynomial operators for parity instead of materialized magnitude [C, S]

Let the target phase be a classical-mask-weighted bilinear pairing. The adjoint/transposed computation may propagate this classical functional *backward through* the decoder recurrence, avoiding the explicit evaluator or root scan. Formulate a reverse/adjoint automatic-differentiation-style propagation over GF(4096), distinguishing linear (Clifford) operations from quantum × quantum field products (charged). Derive on small field cases and compare nonlinear count against current 4,095-point loop. **Critical:** the arbitrary classical mask has 4,095 values; one cannot assume its 70-dimensional syndrome projection is all that matters.

### D17 — Mask-polynomial remainder / quotient functional [C, S]

`dqi-phase` proved a classical polynomial interpolation that eliminates 70 positions using free-to-score classical computation and syndrome Z phases; expected active points are ~4,024. Study whether the remaining classical mask polynomial can be reduced *modulo the quantum locator* with a cheap low-memory transposed product, yielding the total phase without pointwise inverse at 4,025 positions. Compare the full quantum-modulus reduction cost against ~350k scan Toffolis; charge all quantum-dependent operations and clean reversible scratch. Classical preprocessing can use mask bits, but **not** syndromes or true errors.

### D18 — Sparse root extraction plus classical random-mask lookup [C, S]

Only up to 35 error locations matter, versus 4,095 scan positions. Quantum finding roots and indexing the arbitrary classical mask could be expensive, but is worth a measured **end-to-end** screen. Explore additive FFT root factorization, Frobenius splits, product trees, preprocessed mask tables, QROM for 12-bit addresses and batch root evaluation. A 35-root list is not a win if each quantum-indexed classical lookup costs the whole table. Prove each subroutine's reversible storage and handling of variable-degree root counts.

### D19 — Subfield trace functional with shared denominator [B/C, P]

Current 88T body uses 42T GF(4096) product plus a 46T gated cubic trace/inversion phase. Investigate computing the required trace pairing of numerator/denominator **without first materializing** the product: bilinear form simplification, denominator norm, subfield tower, and masked inverse-trace synthesis. Prior 80T evaluator body requires a different upstream representation; copying it into the compact state with a per-point 102T adapter loses. Search for a kernel acting directly on the actual compact `(A,B')` registers, with a measured scratch/gate trade-off.

### D20 — Vectorized GF(4096) operations across a full block [B, S]

Instead of optimizing a generic 42T multiplier, optimize the **entire multivariate map** for a head/tail block using bilinear circuit synthesis, shared linear combinations, small-block tensor-rank search, SAT/MILP or AlphaTensor-style heuristics. For GF(2^12), all linear basis conversions are free to score; only nonlinear AND/Toffoli multiplication terms and quantum storage count. One or two fewer generic scalar gates rarely changes the answer; sharing products across **dozens of correlated field multiplies** could. Verify exact outputs, inverse, scratch and host basis conventions.

### D21 — Search for better field multiplication/interpolation under actual operand constraints [C, S]

42T GF(4096) multiplication is a *found* implementation, not a proved optimum here. Search bilinear maps with nonstandard bases, basis-sharing across successive multipliers, CNOT-heavy decompositions and conditional multiplication by an affine-small operand set. In characteristic two, squaring and multiplication by constants are linear—exploit wherever an operand is provably classical/constant. Do not optimize an unconstrained truth-table multiplier and assume it will fit the decoder's live scratch. The literature includes AlphaTensor-Quantum finite-field multiplication results; test under this benchmark's 12-bit field and Q² objective.

### D22 — Batch inversions / inverse retention across stages [B, S]

The trimmed decoder saves 36T per eligible inverse by holding the inversion chain when live space permits. Investigate shared inverse products, Montgomery-style batch inversion and deferred reciprocal recovery across a block or adjacent heads. Every quantum multiplicand changes by stage, so batch inversion may add product trees/workspace; benchmark full build-and-erase in the real liveness windows. Avoid charging inverse memory only while active; old values must survive until reverse uncompute.

### D23 — Dynamic-programmed expected failure budget and arithmetic schedule [A/B, S]

Extend D06 to a search over *which* gap patterns to support, number of retained exceptional-history entries, coefficient windows, number of updates, block widths, and phase-kernel choice. Use two distinct axes: circuit cost charged for all controlled gates, and failure probability measured under the true squared-eigenvector weight law. Evaluate candidates using the host's exact Clopper–Pearson multiplier; use independent seeds to estimate risk. Optimize Pareto surfaces in `Q,T,expected phase errors` rather than one scalar gate count at a fixed failure target.

### D24 — Exploit cheap Clifford/classical preprocessing safely [B/C, S]

The score ignores CNOT/CZ/SWAP and classical logic gates as Toffolis, though such gates hurt runtime. Use arbitrary classical-mask linear transformations, polynomial basis changes, and syndrome-coordinate maps if they reduce quantum × quantum multiplications. Evaluate memory/runtime effects so 40M additional CX gates do not make the verification campaign impractical. Do not inspect or condition on the hidden error vector or the unknown quantum syndromes as classical values.

### D25 — Interpolate phase only at a carefully selected testable set of points [C, S]

Different from the fixed 70 skip: a classical choice of interpolating positions may minimize total quantum circuit work in a nonuniform scan (e.g. some points have cheaper precomputed field constants or overlap with root-evaluation constraints). If every point costs 11/23/88T identically, choosing *which* 70 to skip provides little value; the purpose is to test genuinely heterogeneous per-position kernels and shared linear maps. Compare all resulting CNOT/Toffoli/scratch costs, not just the number skipped.

### D26 — Approximate decoding by canonical low-complexity subsets [C, S]

Rather than approximating Euclid by dropping generic quotient gaps blindly, design a decoder exact on a statistically selected set of **syndrome forms** that can be identified and handled by cheap reversible invariants, and wrong on the remainder. The allowed error budget should be allocated to structurally expensive edge cases where the omission *removes gates*. Use many independent law samples and adversarial controls, count full data restoration separately from phase parity. This should not be confused with hash/seed grinding or sampling-only correctness.

### D27 — Evolutionary synthesis of reversible decoder building blocks [B/C, S]

Create a library of small proven-invariant word-level transformations (Schur head, local quotient chart, tail MAC, measured erasure, block parity, in-place inverse). Evolve **whole reversible schedules** with exact component costs and liveness at small N, then transfer the same algebraic identities to N70. Reward reductions in `log(R)+2log(Q)+log(T)` and preserve independent correctness constraints. Use self-play across different constructions, not just local peephole edits. Expect the most valuable changes to be representation and state-lifetime changes, not one-gate tweaks.

---

## 3. Chemistry: 28 research directions

**Current anchor:** `chem-census/candidates/structured`, `lambda_eff≈48.56733747`, `T=31,608`, `Q=454`, cost **696.943M**, physically equivalent to the verified single-block R1/B719/C350 8-bit dictionary. It uses J4 / 19 oracles and a nonlinear four-row lookup. The previous J7 / 11-oracle architecture has lower T but greater persistent width. The inner PREPARE cost is **17,333T** and its temporary workspace sets the width peak at **337Q**. With fixed physics, lowering only `lambda_eff` is limited to ~3.7% from this anchor by the proven 46.786 floor.

**Important rank fact:** A legal DFTHC two-body tensor corresponds to a sum of PSD block covariances with total rank <= `R*C` (not the number of dictionary rows). The fresh target's spectral tail requires `R*C>=290` to reach <=20mHa. A smaller R1/B719/C256 shape cannot work at any dictionary; rank-reduced approaches must respect this bound. Realistic target families must also meet quantized vector accuracy, one-body error, normalization/gap, and the exact host QROM widths.

### C01 — Joint nonlinear compression of the J4 four-row bank [A, E+P, leading follow-up]

`chem-qrom` uses 79 retained bank bits and recovers five omitted bits through eight shared products of affine linear forms; it fits exactly in 128Q. Its current expression is `missing(z)=A(z) XOR sum_k output_k*(L_k(z) AND R_k(z))` under 995 row/bit constraints. Optimize **which** columns are omitted, the output basis, `L_k/R_k`, product-output vectors, and row grouping **jointly**; previously the search separated some of these choices. Try randomized factorization over GF(2), SAT/SMT, ILP, rank-revealing greedy removal and local improvement; verify every row exactly. The target is fewer products or more omitted columns *without increasing peak Q*. A one-product reduction is not automatically a task-score gain; compile a whole family and recompute persistent width.

### C02 — Higher-order four-row decoder with controlled cubic/quartic terms [B, P]

The present decoder is quadratic in retained bank bits. Some source bitplanes may compress more strongly using a small number of *cubic* terms or shared nonlinear precomputed monomials. Investigate depth-two AND networks whose intermediate products are recycled via measurement erase. Account for scratch, inverse, and phase repair on arbitrary incoming targets; output parity maps are free Cliffords. Start with one oracle and count actual T/Q before updating all 19.

### C03 — Five-/six-/eight-row lookup with learned multi-output compression [B, S]

The best J is not obviously four. A 2- or 3-bit low address selects 4/8 rows; bank words grow, but compressed feature rank may grow sublinearly if the **data representation is reorganized**. Explore J=5,6,8,10 with nonlinear bank decoding, staged selection and streaming scratch. A naive J7 ANF loses, J13 direct ANF has higher Q, and ordinary J4 two-row lookup loses; avoid those exact implementations. Optimize the *full* `lambda*T*Q`, not Toffolis per oracle. The per-oracle 128Q cap and extra persistent words are critical.

### C04 — Joint low-address grouping, dictionary permutation, and decoder synthesis [B, E+P]

The current J4 compiler gained 34 Toffolis by grouping rows according to three shared linear features, and preserves the tensor by applying the same permutation to U/W. Generic Hamming/lexicographic order and whole-payload repetition cannot lower the J7 count because all relevant payloads are distinct. Search row groups of four/eight chosen to minimize **nonlinear feature rank** and exact decoder complexity rather than Hamming distance. Keep one-body addresses in the correct sector. Search grouping and decoder simultaneously; compile and test every proposed arrangement.

### C05 — Joint basis change across word bits / GF(2) low-rank bank code [B, S]

An invertible output XOR basis is Clifford-inexpensive and may make a nonlinear bank smaller even when it cannot change the **linear** J7 rank49. Use global or per-oracle invertible GF(2) transforms, and optimize the nonlinear decoded bits under the transformed basis. Allow address-conditioned *linear* decoding if its ancillas fit the 128Q cap. Do not claim a random GF(2) basis reduces minimum rank: it preserves linear rank; gains must come from higher-order factoring, zero patterns or grouping.

### C06 — Care-set completion + ANF streamed phase correction [A, E]

The measured care-prefix completion preserves 795 valid rows while choosing behavior at invalid padded addresses to zero high-order ANF terms; seven T per J7 oracle saved. Investigate **jointly selecting invalid-address completions and bank encoding** for J4/J5 so more nonlinear monomials disappear. Make sure the behavior on *all valid padded addresses* remains correct. Reuse the working streaming ANF eraser: it replaces resident monomial scratch by a few live bits, at equal T, allowing tighter 128Q designs. Full forward/reverse random-target checks are essential.

### C07 — Variable ANF bank partition / measurement schedule [A/B, E+P]

Search different high/low address splits and phase-repair orders jointly with compressed bank size. Balance Toffolis for monomials against **public persistent output bits** and peak clean scratch; deleting an output word early may require replay and negate savings. The new report shows J4 is feasible at exactly 128Q; even freeing 1–4 ancillas may unlock a larger nonlinear decoder. Benchmark exact source data, not random tables alone.

### C08 — Construct a single circuit family with shared preprocessing across successive oracle calls [B/C, S]

The interface is cumulative XOR differences: oracle k maps the prior batch words into the next. Could a common address-dependent nonlinear feature register be **persisted** across calls, reused, then erased at the end? The verifier calls separate circuits and the host charges persistent registers; any scratch required zero after *each* oracle may prohibit persistence. First read the precise per-oracle contract. Only explore cross-oracle caching if the host permits that state and it improves the charged total. Otherwise exploit **compile-time** common subexpressions across oracle definitions, not fake live caches.

### C09 — Special one-body sector compression using its independent tolerance [B, S]

The accepted error1 is ~11.995mHa, leaving ~8mHa before its separate 20mHa cap. The 76 one-body ANGLE rows may admit cheaper quantized words or shared code patterns by allowing small physical shifts, while keeping error2 fixed. Uniform one-bit rounding failed disastrously (~123mHa error1); try per-vector sensitivity, eigenvalue-aware matching, sparse candidate moves and constraints on the actual one-body matrix. Score whether this reduces *compiled* lookup T/Q, not merely angle entropy. If no decoder compression results, stop; the spare one-body error cannot subsidize error2.

### C10 — Physics/compiler co-design: optimize dictionary *and* Boolean encoding complexity [B/C, radical]

Optimize quantized row words to make J4 bank features lower-rank while constraining the true two-body residual and `lambda_eff`. Example: introduce a penalty on the number of nonlinear decoder products, or on projected bank-bit rank in batches, not generic Hamming distance. Alternate small word moves, W/H2_I refit, row permutation and exact circuit synthesis. Use the accepted physical basin initially to prevent losing feasibility; require every candidate to improve *predicted actual* `lambda*T*Q` and pass final native. This is a different objective from pure tensor fitting.

### C11 — Quantization-robust dictionaries with explicitly conditioned feature spans [B/C, central radical]

`chem-tree` found B351 continuous factor blocks with condition numbers ~1525/2210 and weight Frobenius norm ~55, versus incumbent ~172 and ~3.72. Eight-bit rounding destroys a few nearly singular directions: worst principal angles between continuous and quantized spans reach 82–88° even though median span cosine is near 1. **Optimize conditioning directly**: penalize inverse smallest singular values of the *weighted projector feature* Gram, leverage concentration, and the spectral norm of quantization-induced perturbations. Impose usable normalization as a constraint *during* dictionary construction. Then quantize and refit. The aim is a *different well-conditioned basin*, not to round the fragile C189 basin harder.

### C12 — Train with actual integer Givens decoding and structured perturbations [B, P]

Continuous fits of R2/B351/C189 pass (19.979mHa) but rounded words have >140mHa union-span error. Gaussian noise/STE alone has been unsuccessful in prior rounds. Develop a **quantization-in-the-loop optimizer** using exact decoded words in forward passes, explicit candidate accept/reject, antithetic half-LSB perturbations as sensitivity estimates, and periodic integer-word reinitialization. Train the actual physical tensor with a spectral conditioning penalty and low Lambda. Test on the *accepted R1/B719* as a control: the method must preserve that basin before concluding it improves harder shapes. Do not use one failed STE run as a theorem about QAT.

### C13 — Residual-guided multi-atom exchanges with annealing and low-Lambda constraints [B, P]

`chem-discrete` tested correlated word moves/atom replacements and found that matched spectral refits usually beat complicated moves on its specific B479 start. Extend with **nonlocal supports**, pair swaps across ensembles, residual-mode quantized atoms, and annealed acceptance on *both* error and Lambda. Crucial: a candidate that reaches ~20mHa with Lambda 100 is unlikely to win. Store the actual Pareto frontier `(error2, Lambda, T, Q)` rather than a single best error. Give each move family enough CPU/GPU work to compare at matched eigensolve counts; abandon only the tested family, not integer combinatorial optimization.

### C14 — Low-Lambda pruning along known feasible physical generators [A/B, P]

The banked incumbent at B719 is extremely close to the physical threshold. `chem-census` removed 16 rows to B703 and repaired actual words to 19.982mHa, but Lambda deteriorated so its PASS cost was ~807.5M. Instead of free weight refitting after pruning, carry the **original 350 physical generator matrices** as optimization constraints, or use covariance/column gauges preserving their low-norm structure while atoms are removed. Select deletions by *marginal physical error and projected normalization increase*, not leverage alone. Search swaps/new atoms near each width cliff and a constrained homotopy path toward B639/B511.

### C15 — Integer-native dictionary construction (growth) with covariance-matched atoms [B, P]

`chem-census` grew quantized residual-mode dictionaries and was decisively worse than reinserting incumbent atoms for its tested pools. Consider a new atom generator that **matches the tangent or projected Hessian of the low-Lambda covariance**, not simply the leading residual eigenvector. Use multiple candidate directions jointly, QR/Gram-aware orthogonalization, target signed tensor including its small negative modes, exact word quantization and a constrained rank-C PSD solve. Compare against sequential reinsertions and match resource budgets. Treat the saved B511/B671 growth failures as local capacity bounds, not B-wide bounds.

### C16 — Optimize the B/C/R shape for host QROM width discontinuities [A/B, E]

The host selects internal QROAM exponent `k2` automatically; it is **not** a participant-controlled knob. Thresholds are not always monotone in B or brot: B513->512->511 changes Q478->462->460; R4/B180 actually has k2=4 for both 7/8 bits, contrary to a previously optimistic example. Create an accurate shape × bit precision × J × actual-oracle compiler frontier (using measured or defensible bounds) and overlay feasible quantized dictionary trajectories. Try even counterintuitive *padding* to cross a favorable `k2` cliff if allowed. Do not spend large fitting budgets on a shape with little or negative conditional score upside.

### C17 — Minimize C toward the universal rank threshold of 290 [B, E]

At R1/B719, C320 is excluded **for the fixed accepted dictionary** (~21.36mHa), and the global necessary rank excludes `C<290` at R1. A new dictionary might approach C290–310. Optimize dictionary geometry and column covariance together, with careful conditioning near the spectral tail. Calculate host T savings and Lambda trade-offs. Do not repeat the provably impossible standalone C224/C256 plan or confuse `B` vector count with `C` factor rank. A simultaneous B and C reduction may cross multiple resource discontinuities.

### C18 — Block-PSD relaxations as *design guides*, not just kill certificates [B, E+P]

`chem-psd` and `chem-discrete` built accurate primal/dual PSD certificates for fixed dictionaries; these are powerful tools. Use the **dual residual eigenmodes** to generate new quantized projector vectors and guide which ensemble should receive them. For each fixed dictionary, separate union-span, block-linear, full-PSD, rank-C and exact decoded errors. Optimize *the gap between these relaxations* rather than blindly running L-BFGS. A dual bound >20mHa justifies changing the **words/dictionary**, not abandoning its R/B/C shape.

### C19 — Learn block assignments / overlapping subspaces instead of joint diagonalization [B/C, P]

Pure JAD and R4/B180 substantially failed continuous accuracy, even after long refits; the matrices do not share enough eigenvectors. Revisit multi-block factoring with **nonorthogonal overlapping ensembles**, shared high-leverage atoms, unequal physical error budgets, learned mode-column partitioning, and constrained covariance solves. Do not require a nearly common eigenbasis, and do not confuse the discredited *orbital-disjoint* split (which omits cross-orbital correlations) with legal full-orbital overlapping dictionaries. Use resource-first shape selection and a PSD-optimality screen before running expensive quantization.

### C20 — New positive-cone decomposition of the signed target with sparse rank-one atoms [C, S]

The verifier's two-body tensor is approximately PSD in symmetric pair space, with a small negative spectral component (~0.298mHa). Consider a **conic decomposition**: a collection of rank-one *orbital* projectors and PSD covariance modes that approximates its positive part in a quantization-robust fashion. Use column generation: price new dictionary vectors against the current convex dual tensor, insert a bundle, refit covariance, prune. This is analogous to cutting-plane/Frank–Wolfe atomic methods but over *rank-one orbital projectors*, not arbitrary symmetric-pair eigenvectors. Keep normalization penalties from the beginning. The convexity holds for weights at fixed dictionary, not jointly over vector directions.

### C21 — Joint LCU normalization and discrete geometry with an active inequality barrier [B, E+P]

The universal Lambda bound is only 3–4% below current; that **does not mean** smaller shape candidates can tolerate Lambda doubling. Solve `min lambda_eff*T*Q` subject to exact residual<=0.020, strict gap, and legal codes; use a differentiable augmented Lagrangian or a feasible-only filter/trust method, with exact host metrics at checkpoints. Penalties on the **reported Lambda** are distinct from arbitrary L2 norm penalties. Preserve feasible candidates even when a numerical optimizer explores invalid intermediate states. Consider multiobjective Pareto evolution when the feasible set is thin.

### C22 — Orthogonal rotations of factor columns / H2 offsets [A/B, E]

Column rotations preserve the represented generator covariance when weights and identity offsets transform together, but change the L1-type normalization and lookup control opportunities. Rotate within blocks to concentrate coefficients, produce simple signs, or reduce cancellation while tracking `E_gap`. Prior sign-mixing alone improved only ~0.001% beyond control; use this **as an enabling gauge after a dictionary/shape change**, not as a standalone road to 2×. Verify the exact augmented Gram and one-body offset contraction.

### C23 — Per-address Givens representability with a *learned shared tree* [B/C, E+P]

`chem-tree` tested 17 distinct binary amplitude lattices and ~900 ordered trees; tiny gains, major errors remained. But it also found quantization-sensitive directions and mathematically equivalent tree signatures. Start from **new well-conditioned dictionaries trained inside an alternative binary lattice**, not a fragile continuous dictionary encoded under a new tree after the fact. Jointly optimize tree, integer words, weights and normalization, and compare at matched complete compute budgets. Tree mutation alone on the old exact words is deprioritized.

### C24 — Legal angle coding / signed gauges / nine-bit trade-offs [B/C, E]

The host semantics enforce `theta=pi*word/2^(brot-1)` for uniform `brot`, with one shared valid rooted tree. Eight-bit sign/gauge equivalence does **not** expand the representable projector set: the `chem-tree` report checked this. A 9-bit shifted grid has more physical resolution but can trigger k2=5 and near-double full width (e.g. Q778 for R2/B351). Search **shapes** for which the nine-bit precision helps and the preparation k2 threshold remains favorable; alternatively exploit special eight-bit subsets to compress the lookup, not to claim improved angle resolution. Charge any decoding and bitword width honestly.

### C25 — Address-aware compression of one-body eigendata separate from two-body rows [B, S]

Currently 76 one-body angle rows precede the `R*B` two-body rows. If one-body vectors share structured relations, build a separate affine/low-degree representation whose values are integrated into the same mandatory oracle family. One-body errors are checked independently and eigenvectors are matched by the host; ensure reordering does not alter their semantics. A cheap one-body block could free scratch or reduce total Toffolis even when the two-body table is dense. Test on the 76 actual rows before inventing constraints on 795 rows.

### C26 — Physically structured dictionaries that *also* have compressible angle words [C, S]

Design dictionary atoms as products of sparse orbital cluster amplitudes, near-Hadamard/simple angle states, or a mixture of generic and low-complexity atoms, with **all 76 orbitals allowed**. Then fit PSD weights and measure physical cost. Pure locality/disjoint blocks failed due to missing Coulomb correlations; use sparse atoms only as part of an **overcomplete global** dictionary. The possible payoff is a cheaper QROM (many zero/simple word bits) and robust quantization, but the risk is severe loss of representational capacity.

### C27 — Alternative physical factorization search from accepted generator columns [B, P]

The incumbent gives 350 symmetric generator matrices. Rotate those generators orthogonally (which preserves the unquantized tensor when offsets rotate correctly), cluster them by **joint physical covariance / leverage**, and rebuild a compact dictionary that need not be a direct partition of old rows. JAD's commutator diagnostic was poorly predictive; replace it with block-PSD approximation quality and conditioning. Use spectral line search plus exact quantized word generation, not merely smooth orthogonal diagonalization. Try learned overlapping ensembles and rank-C fitting.

### C28 — End-to-end compiler-and-factorization evolutionary search [C, S]

Evolve tuples `(R,B,C,brot,J,tree,quantized_U,W,H2_I,oracle architecture,order)` with a hierarchical evaluator: fast rank/span and resource feasibility, then accurate PSD/physical metrics, then native oracle compilation, then complete verification. This searches **representations** instead of local weights. Use a genetic/evolutionary archive preserving diverse Pareto alternatives and periodic reseeding from successful physics or compiler candidates. A GPU optimizer that only minimizes residual can repeatedly produce noncompetitive high-Lambda models; the expensive objective must influence candidate selection from the start.

---

## 4. Inverse square root: 22 research directions

**Current anchor:** 90Q, 1,782.00009 executed Toffolis, 34 synthesized rotations, ~250.85 equivalent synthesis work, cost **182,956.52716**. The 90Q peak occurs in the **first radix-four row of the final Horner product**; normalized input, derivative/slope, late-loaded constant, curvature and carry storage are already aggressively aliased. The chosen 269-leaf prefix table is independently rediscovered by `inv-dyadic`. Any win must alter nonlinear arithmetic or its scratch/phase lifetime, not just call another interval partition "adaptive".

### I01 — **Build the untested narrow fused-square quadratic** [A, P, top inverse lead]

`inv-census` found an exhaustive-numerically passing pipeline that keeps the linear term exact but rounds the local coordinate only **inside the quadratic term**. With signed `z`, set `u=floor((z+16)/32)` and `h=round(u²/64)` using the report's precise integer rule, then

```
num = c0*2**23 - c1*z*2**11 + c2*h*2**16
phase_units = num * 2**(s+4)  (mod 2**64)
```

Its exact 3,133,452-input model passes with worst error/tolerance ~**0.804696** and absolute margin ~7.48e-12 turns. The rounded coordinate needs **nine unsigned magnitude bits including endpoint 256 plus an extra sign/guard**; the square is ~11 bits. **First test:** implement a *specialized reversible square* and fused `c2*h` phase/accumulation with explicit live-register reuse. Compare at 78–90Q; exact output model + small native controls before a full run. Do not round the linear term or forget endpoint/sign behavior.

### I02 — Fused exact-rational quadratic with high synthesis margin [B, P]

`inv-phase` proved the exact-rational quadratic representation

```
P_exact = C0 - C1*z/2**12 + C2*z*z/2**23
phase_units = (C0*2**23 - C1*z*2**11 + C2*z*z) * 2**(s+4)
```

passes all inputs at worst relative-error ratio ~0.62666 and minimum synthesis margin **1.6008e-11 turns**—~14.8× more margin than the incumbent. This is *not* a built competitive circuit; straightforward term-by-term phasing adds too many rotations. First test: synthesize the **numerator's required phase bits directly** with a narrow fused square/middle product and a **single shared exponent shift**, not 201 independently exponent-routed rotations.

### I03 — Specialized squarer with linear diagonal terms [B, S]

For binary `z`, use `z² = sum_i z_i*2^(2i) + 2*sum_{i<j} z_i*z_j*2^(i+j)`; diagonal products are linear because `z_i²=z_i`. The best reversible squarer may need fewer Toffolis than a generic variable×variable multiplier, especially when high/low output columns are provably discarded after rounding. Derive a truncated square bit-heap and use measured temporary-AND erasure. Optimize bit truncation jointly with exact output-phase tolerance; validate signed two's-complement, guard, and carry semantics.

### I04 — Transposed/middle-product quadratic correction [B, S]

Only a subset of the final product bits contribute to rotations, and their carry effects can sometimes be reconstructed from a shorter middle product. Build a specialized `c2 * h` or `C2*z²` accumulator for the necessary output columns, with cheap linear maps and low-ancilla carry handling. The existing Horner pipeline already uses radix-four Booth and unequal-width adders; the new improvement must be a *different partial-product dependency graph*, not generic schoolbook multiplication. Measure full compute/phase/uncompute, not only forward arithmetic.

### I05 — Normalize exponent information in the input register [B, P]

`inv-census` constructed an injective map over all original inputs: `E=M` for k=0 and `E=M-2^21+2^(k-1)` otherwise; low marker bits permit recovering k by a trailing-zero count. A complete circuit was not built and the coordinate changes. Derive a reversible normalization/denormalization that reuses these marker bits instead of reserving five exponent flags at the Horner peak; reconstruct k on demand. Price the trailing-zero comparator and any local-coordinate repair. A five-qubit saving is meaningful at 90Q only if it costs <~100 Toffoli-equivalent extra work.

### I06 — Phase exponent routing shared across arithmetic, not repeated across rows [B, P]

The direct-row phase implementation saved last-product arithmetic but spent **201 rotations** and hundreds of exponent shifts, scoring 318k. A Walsh selector needed ~307 rotations, and `inv-census` proved 10-character exact Walsh support cannot cover the tested per-bit multiplexer. Revisit only if a **single common phase bank** or different joint bit encoding lets *all* rows share the exponent shift. Do not rerun independent rotations on each row or individual per-bit Walsh interpolation unchanged.

### I07 — Signed phase-word wrap corrections with near-zero Toffolis [B, P]

`inv-phase` derived exact sign handling by CNOT-complementing a phase bank and using one correcting rotation; the approach is valid but total rotations dominate. Explore transforming the signed output word into a globally positive biased phase word that shares the same shift across rows, or fold wrap correction into a small number of existing rotations. Build exact integer proofs for all exponent s before compaction. Count synthesis work from the whole circuit, not from the number of arithmetic Toffolis saved.

### I08 — Phase-oriented Booth bit heap instead of final Horner accumulator [B, S]

The incumbent's 1,046 static forward+reverse Horner Toffolis are central. Try a **carry-save/bit-heap** representation that keeps a pair of sparse bitplanes or signed digits through both polynomial stages and collapses only the phase-relevant columns. The bit heap may use more qubits than a ripple accumulator; optimize `Q*T` with a realistic liveness profile and measurement-erased carries. Reject if the second independent plane increases Q more than gate savings repay. Compare to the 90Q incumbent's exact scratch overlap.

### I09 — Direct approximate phase polynomial in coarse+residual form [C, S]

Approximate `1/sqrt(x)` as a coarse prefix-dependent phase plus a residual generated from a narrow local coordinate, with both designed to use the **same** small set of rotations. The prior pure Walsh/Boolean ANF representations were uncompetitive because they needed too many phase terms; a new approach must exploit monotonicity/curvature and precomputed **shared** high bits, not generic Fourier sparsity that is absent. Certify every x against the host integer tolerance and all synthesis margins.

### I10 — Piecewise rational/minimax approximation optimized for *reversible* arithmetic [C, S]

Investigate Padé or rational functions on normalized intervals only if numerator/denominator evaluation plus inversion can be cheaper than the current quadratic Horner, or if a denominator takes a specially cheap form (power of two, small correction, fixed affine term). A generic Newton/Goldschmidt reciprocal was already tried and lost. Prototype integer approximants under fixed-point truncation, then prove the denominator's reversible inverse is cheap before developing a complete circuit. Treat superior continuous approximation as insufficient evidence.

### I11 — Early phase plus measurement-erased carries with correct coherence [B, P]

`inv-width` built segmented addition with a retained carry and HMR phase repair, saving some width but increasing Toffolis. The last Horner multiplication has a **multi-bit signed carry**, not a single one: after a split at low bit 14 it needs ~6 signed carry bits. Search a special representation in which residual carry information is retained compactly or reconstructed during inverse execution. Check that measuring it doesn't leak x-dependent information into an uncorrected phase. At each split report both *actual retained bits* and total gates.

### I12 — Reversible low-ancilla adder adapted to dirty borrowed bits [B, P]

Existing dirty-borrowed adders reduce qubits but doubled or tripled gates in tested independent pipelines. Instead specialize the existing **short uneven-width** adder to its exact operand promise: e.g. known-zero high bits, fixed sign, known carry ranges, and distinct polynomial row widths. Search SAT/peephole or local scheduling for those promises; control arbitrary dirty-borrow values. A generic Cuccaro/UMA swap is not enough, and the incumbent already uses input/curvature/output wires as borrowed scratch.

### I13 — Joint search over precision, coding and actual synthesis margin [B, E]

The inverse task isn't scored solely on maximum phase-error ratio; `min_x(tolerance(x)-error(x))` sets the cost of **all** remaining rotations. Some more-accurate polynomials reduce synthesis cost, but extra table leaves or coefficient bits may outweigh it. Optimize exact integer coefficient fitting against **true end-to-end circuit score**, not an RMS or uniform max-error target. Search asymmetric margin allocation across normalized intervals and output bits, and test whether 1–3 coefficient changes can enlarge the *minimum absolute* margin without adding Toffolis. Use the existing 269-leaf source as a calibration and don't expect huge savings from synthesis alone (~251/2033 work units).

### I14 — Input-domain reachability and parity-specific circuits [B, E]

The normalized mantissa's legal values differ by exponent parity: some values are never reached. The existing integer fitter exploits these restrictions; it is worth testing more aggressive **region-specific** fixed-point transformations that are valid only on reachable subspaces and are coherently controlled by quantum exponent flags. Prove every x in the finite original domain and test branch cleanup. Never set a register by assuming the actual quantum x is a classical constant; only algebraic reachability identities can justify omitted gates.

### I15 — Address table factored into affine base plus sparse curvature residual [B, E]

Round 4 already removed a cheap affine function of interval address to reduce stored constants. Search a jointly fitted **low-degree address predictor** (Clifford-derived bit parities, affine or small nonlinear terms) plus a much narrower QROM residual. Measure both the cost to construct predictor and the reduced bank; avoid replacing a 267T coefficient lookup with a 500T arithmetic predictor. The inverse table's apparent 512 rows already compresses to 269 leaves, so row count is not the correct baseline.

### I16 — Use alternate normalization bases (radix 2/4/8, exponent/parity folding) [B/C, S]

The current normalization is very inexpensive (~210T, five flags), but a different radix changes coefficient value ranges, necessary output shifts, and the Horner coordinate width. Explore radix-eight/octave or mantissa mapping into `[1,2)` / `[2,4)` with exponent absorbed into a small fixed phase code. The goal is to improve the **downstream multiplication**, not to save 20T in normalization while doubling lookup rows. Determine if a new base removes the 10-bit output headroom without introducing an expensive exponent mux.

### I17 — Small-domain exact lookup for singular input regions [C, S]

The smallest x values may force a disproportionate numerical precision or absolute synthesis-margin constraint. Consider carving out the first few dyadic regions with exact direct phase constants or bespoke low-rotation mini-circuits, allowing the main polynomial to run at narrower precision elsewhere. Quantum selection/gating still costs gates; use the complete input distribution (invsqrt is **exhaustive**, so no rare-case probabilistic dismissal). Verify whether the global minimum synthesis margin actually occurs in that carved region before spending gates.

### I18 — Variable polynomial degree with shared reversible schedule [C, P]

The `inv-dyadic` agent tested mixed linear/quadratic tails, losing because extra leaves and a shared quadratic schedule remained. A potentially different approach is to select **different complete arithmetic paths** for linear and quadratic regions with a cheap reversible mux and *shared* registers; charge both controlled circuits, not their input occurrence rate. Design only if most of the normalized domain can omit a substantial nonlinear stage, otherwise the fixed overhead dominates. Search physically meaningful prefix regions and proper absolute tolerance.

### I19 — Constant-coefficient multiplication with addition chains or specialized NAF [B, S]

Although coefficients are loaded from a table, their small bit widths and possible common bit patterns may permit cheaper multiplications than generic Booth Horner. Search signed-digit/NAF decomposition after a different angle/base encoding, reuse shifted input partial products across coefficient families, or precompute small common products of the local coordinate. In the quantum circuit coefficients remain *quantum-dependent on x*; treating them as compile-time constants is invalid. Test an actual shared reversible arithmetic gadget at all table addresses.

### I20 — Shift/rotation synthesis from bitwise phase changes without materializing a full word [C, P]

The full result after scaling is a phase, not an output. Investigate whether a fused multiply-add can emit one or two groups of phase rotations while freeing dead partial columns, without paying the 201-rotation penalty of direct row phasing. The condition is shared exponent routing and few total rotations; a proof of integer phase additivity is useful, but every independent output rotation is charged by host synthesis. Compare after native margin recomputation.

### I21 — Global reversible circuit superoptimization of the actual incumbent bottleneck [B/C, S]

Use a symbolic/QSAT/ABC-guided search on the first radix-four final Horner row or its unequal-width adder. Borrow exact semantic promises on operand sign, masks, carry range and scratch cleanliness; search maps that use 3–10 fewer qubits or 100–300 fewer gates. Small changes can improve overall task cost by ~5–15%, unlike blindly changing approximation formulae. Prove circuit equivalence with independent integer models and native random dirty-input checks. Do not feed all 69M operations into a generic optimizer; optimize a constrained subroutine and reinsert it.

### I22 — Mathematical representation chosen for optimum quantum *space-time*, not classical FLOPs [C, S]

Search beyond polynomial interpolation: algebraic differential updates, Chebyshev/Clenshaw with truncated phase accumulation, a short rationalized Padé correction, fixed-point bipartite phase tables, or a bit-sliced reciprocal-sqrt digit recurrence. For each, create a **tiny symbolic cost model first**: coefficient load, nonlinear products, phase bank, rotations, scratch, reversal. Reject families only when a credible optimistic lower-bound screen already exceeds 182,956; otherwise implement the cheapest primitive rather than a whole circuit immediately. This is where a score-20 agent should look if conventional Horner is truly saturated.

---

## 5. Cross-cutting methods that may unlock several tasks

These are not separate physics tasks; they can make a radical idea *actually buildable*.

**X01 — Quantum circuit superoptimizer guided by actual live-register intervals [A/B].** Design a resource IR that annotates operations with data dependencies, clean/dirty ancilla promises, measurement erasure, and Clifford-only basis maps. Enumerate schedules and register aliasing, **including** reverse cleanup and phase correction; optimize a peak-Q-aware objective. The existing DQI and inverse emitters perform local compaction but may not be globally optimal under alternative algorithms. Test on a known circuit to quantify the IR's own accuracy before trusting projections.

**X02 — Charge nonlinear rank, not classical FLOPs [B].** GF(2) basis transforms, affine parity networks, CNOT/CZ and classical-mask processing may be very long but contribute little or no Toffoli cost under the present objective. Design arithmetic circuits directly for **multiplicative complexity** and available clean ancillas rather than importing optimized CPU/GPU algebra. Also measure simulator size/runtime so the circuit is still verifiable overnight.

**X03 — Automated construction of reversibility proofs [B/C].** Given a proposed destructive update, search for a recoverable invariant or an explicitly stored minimal history tag. Prove the resulting update injective/permutative on the promised input domain and test inverse-first composition. This is especially important for DQI zero pivots and quantum polynomial quotients. Branching based on a quantum value requires coherent gates; it cannot be treated like a CPU if-statement.

**X04 — Error-budget and numerical-margin co-design [B].** DQI can tolerate statistically priced failures; inverse cannot tolerate any checked error and pays global synthesis-margin costs; chemistry has two independent 20mHa residuals and strict SOS floors. The agents should use **three different notions** of approximation. Do not import one task's trade-offs into another without verifying the contract.

**X05 — Structural counterexample mining [A].** Use independent classical oracles and native kickmix components to identify the **first divergence** from the intended mathematical invariant. Minimize failing instances; add them as permanent regression tests. For optimizer failures, compute per-shape relaxation bounds, conditioning, and exact quantized physical residual to distinguish numerical stagnation from capacity. The round-5 DQI repair succeeded because it diagnosed *why* the 380 errors occurred.

**X06 — Hierarchical physical/circuit evaluation [A].** First mathematical identity, second independent classical simulation, third native small components, fourth exact emitted Q/T, fifth full native verifier, sixth packaged metric. This allows long exploratory research without treating projections as accepted scores. When an accepted candidate is found, immediately freeze bytes, independent-salt verify and write a machine-readable selection record.

**X07 — Port only the useful part of external research [A/B].** The public `andigu/uncompute` RS decoder is much cheaper in T but too wide under this host Q² objective. Its algorithms, normalized field maps, history reconstruction, and masks are useful **components**. Likewise, published Coulomb/HHL approximation and Fe2S2 lookup compilers have different tolerances/oracle shapes. Cite the source/license and adapt the algebra, not the score projection.

**X08 — Public-notebook mining with source/implementation verification [A/B].** Kaggle's public code listings currently surface titles associated with direct polynomial phasing, compact arithmetic, and banked QROM; those titles **do not establish leaderboard score or implementation quality**. If public notebook sources are accessible, inspect reproducible code (with attribution and licensing), identify techniques not represented in the local backlog, and re-evaluate under the *exact* host verifier. Don't assume the leaderboard leader's private method is public. Prioritize concrete implementations over promotional descriptions.

**X09 — Distinct representation populations to escape LLM local optima [B].** Maintain simultaneous candidate families that share only the native verifier, not their mathematical assumptions: e.g. DQI C-fraction vs Wronskian vs mask functional; chemistry one-block quantum-code-friendly vs robust multi-block; inverse fused-square vs rational-approx vs alternate phase bank. Archive diversity prevents agents from mutating the current best circuit until every branch looks alike. Keep explicit negative results and attempt one structurally distant reset when the Pareto frontier stagnates.

**X10 — Search for resource-model discontinuities and cross-task scoring arbitrage [A/B].** The host cost model changes at QROM block widths, batch sizes and accepted repetition counts. Sweep all **legal** shape/control variables and flag discontinuities. Even an apparently counterintuitive larger B, an additional classical bit, or fewer loaded words per oracle can lower total `T*Q`. This is ordinary model optimization, **not** modifying the scoring system or inputs.

---

## 6. Salvage ledger: don't accidentally re-run discredited implementations

| Tested approach | Credible negative finding | What could still salvage the mathematical idea |
|---|---|---|
| DQI standalone half-degree Psi | 1,383Q implementation ~1.71× worse than original 961Q reference even with free evaluator conversion | Produce a *different decoder state* that evolves Psi/H directly, or alias its workbank safely into consumed syndrome words |
| DQI in-place F postprocessor | Measured fast modular MAC lifetime 93,408T and >2,900Q component | Derive a genuinely in-place reciprocal/phase-aware recurrence; never allocate literal old-F/new-F banks |
| DQI generic C-fraction serial decoder | ~1.123M T core build/erase, unsupported zero-head charts too common | Blocked coefficient extraction with cheap head/tail update and coherent exceptional-chart repair |
| DQI direct Wronskian/H sandwich | 46T point body and H preparation work, but chart errors and integration costs erase margin | Update `(B,H)` during decoding, cheap quotient-driven transfer, share phases with locator |
| DQI old BM template reuse | 840 input + 432 separate locator bits -> at least 1,272Q even without other ancillas | New representation that stores locator information *in* the consumed syndrome bank, not shortening existing lifetime slightly |
| DQI separate carry-chart dispatch | Multiple charged moment maps cost +326kT *before* full correction | One shared nonlinear moment map with Clifford routing; learned factored transfer representation |
| DQI matrix/Bennett transfer | Explicit 2×2 matrix coefficient bridge costs 4,356T per generation and >1,107Q with other state | Avoid materialized matrix entries entirely: factorized quotient transfer or shared Schur head |
| DQI fully correct gap decoder alone | Verified 1.113× improvement, core 824,740T / 1004Q | Partial-correction resource/failure optimization; cost-aware algebraic simplification; combine with trim where compatible |
| Chemistry trivial rank pruning | C<290 impossible for **any** legal PSD factorization in this host; old dictionary at C320 fails | Jointly change vector dictionary and C in range >=290; reduce B, QROM or J as well |
| Chemistry continuous R2/B351/C189 then direct rounding | Continuous 19.98mHa; rounded dictionary's union span >140mHa | New integer-native **well-conditioned** dictionary, not stronger nearest rounding of this basin |
| Chemistry fixed B479 saved dictionary | Convex block-PSD lower bound 20.156mHa (>20) | Change the exact integer word table or geometry; don't just refit W or increase C |
| Chemistry pure JAD R4/B180 | Full-rank PSD lower bound ~50.5mHa on tested dictionary | Overlapping, nonorthogonal dictionaries and different generator geometry; don't insist on shared eigenvectors |
| Chemistry global normalization-only search | Universal `lambda_eff >46.786` prevents >~3.86% from current physical normalization | Reduce T/Q; exploit normalization as a constraint during compression, not standalone breakthrough |
| Chemistry simple repeated-subtree pruning/row reorder | All 795 source full payloads and pair XORs distinct; naive permutations save 0T in old compiler | Higher-order bitwise algebraic factoring, nonlinear bank encoding, J4+/care completion, jointly chosen address groups |
| Chemistry eight-bit signed gauge variations | 128 high-bit gauges had unchanged J7 bank rank49, same count | Change compiler nonlinear basis/physical dictionary or shared tree, not another pass over equivalent high-bit signs |
| Chemistry plain alternative tree after nearest rounding | Several hundred tree proposals mildly improve B351 but leave 54mHa+ errors | Construct a *well-conditioned integer dictionary directly* within a new binary lattice |
| Inverse dyadic prefix repartition | Independently recovered near-identical 269-leaf incumbent; widths/gates tie | Change arithmetic representation (fused-square) or phase routing, not one more interval partition |
| Inverse direct seven-row phase | ~201 rotations, 318k full native cost | One **shared** exponent/phase bank with middle-product grouping; no independent exponent shift per row |
| Inverse biased split phase | 84Q but 63 rotations / extra table replay, 231k cost | Fold constant/residual into one rotation bank and avoid redundant curvature replay |
| Inverse Walsh bit multiplexer | >300 rotations; exact lower bound for **that** 10-character bitwise representation | Change joint bit encoding or approximate phase in a different basis; don't repeat 10-character search |
| Inverse dirty-borrowed generic adder | Quibits drop but nonlinear gate count spikes | Specialize to exact truncated-product carries; fuse with squarer, rather than borrowing dirty bits in generic serial Horner |
| Inverse segmented adder | Retained multiplication carry needs multiple signed bits, not one | Change carry representation/phase schedule or exploit narrow fused-square calculation |

**Interpretation:** This table is a map of *failed implementations*, not a ban on their associated mathematical problems. An agent claiming global impossibility must identify an actual proof whose assumptions cover the new representation. Likewise, an agent resurrecting a failed implementation must state what crucial constraint it has changed.

---

## 7. Autonomous `/goal` research protocol for unattended runs

> **Suggested operational constraint:** Treat this document as an idea portfolio and a source of hypotheses, not a mandatory sequential checklist. The agent has autonomy to select and combine ideas, but should not keep applying the same local optimizer indefinitely or stop after one uncompetitive prototype. The targets 13, 17 and 20 describe aspiration; all claimed scores must come from actual native verification.

### 7.1 How to start and keep state

1. **Read repository instructions and immutable evidence.** Discover root `AGENTS.md`/skill instructions and selected source artifacts. Enumerate completed directories under `initial_experiments/5/`; distinguish *finished* from *concurrently active*. Do not peek at active agents' workspaces unless the user explicitly authorizes collaboration.
2. **Freeze exact task-level baselines.** Copy or reference immutable candidates: inverse round-4, gap-aware DQI, structured banked chemistry. Record SHA256, full verifier output, default and independent-salt result where available. The older assembled archive may not contain these exact best components.
3. **Build a shared-format scoreboard in your *own* workspace:** per candidate `task`, `source`, `hash`, `verified_status`, `Q`, `T_mean`, `R` or `lambda_eff`, `task_cost`, `improvement_from_best`, `combined_score_if_substituted`, and `artifact_paths`.
4. **Inspect the source of the current bottleneck.** DQI: core/scan/flags/phase; chemistry: inner PREPARE/ANGLE/persistent words/normalization/physical residual; inverse: final Horner/phase routing/rotation synthesis/peak live workspace.
5. **Choose at least three structurally different hypotheses** whose mathematical inputs and circuit resource model differ, plus one low-risk component experiment. Do not just pick four parameter values for one algorithm.
6. **Use a source-of-truth evidence ledger** (see 7.5). Classify every result: mathematical identity, independent model, native component, whole-circuit estimate, native full pass, packaged metric. Claims of score improvement require a complete verifier PASS and a lower score for the *same task objective*.

### 7.2 Time-budget and research persistence policy

- An agent in `/goal` mode should continue investigating while it has useful avenues, rather than stop after 15–40 minutes because the **first** implementation lost. The user will be unavailable until tomorrow evening; maintain resumable checkpoints and readable progress.
- **At least once per meaningful research phase**, reconsider premises. Ask: “What assumption forces this design to be expensive? Can I remove that assumption without violating the verifier contract?” Candidates: polynomial state, auxiliary history, circuit direction, normalized coordinate, tree, integer representation, word grouping, scratch schedule.
- When a cost screen rules out one branch, write its precise resource inequality and pivot to another family. A failed *first* branch is a signal to update the portfolio, not evidence that the requested score is unreachable.
- Avoid brute-force GPU/CPU loops with no evidence of decreasing objective. Prefer diverse restarts/representations, cached exact targets, structured matrix operations, and checkpoint histories. A long run is warranted only when the actual objective is progressing or the search is exploring *new territory*.
- Do not invent a time limit or a hard number of experiments as an excuse to stop. Equally, don't burn shared resources on a mathematically excluded shape. When all locally viable branches are exhausted, produce an extensive next-research handoff rather than asserting a global ceiling.

### 7.3 Four independent portfolio briefs (optional specialization)

For four autonomous agents targeting different score levels, specialization reduces correlated local optima:

| Portfolio | Suggested directory under `initial_experiments/5/` | Main search space | Secondary backup |
|---|---|---|---|
| **Alpha: decoder representation** | `goal-dqi-repr` | D01–D04, D12–D18, C-fraction/Wronskian/partial inverse; seek 2–3× DQI | DQI partial-gap Pareto optimization and trim integration |
| **Beta: compiler / physics co-design** | `goal-chem-compiler` | C01–C10; exact nonlinear bank compression and QROM+factorization coupling | small validated improvements to selected chemistry |
| **Gamma: new chemistry geometry** | `goal-chem-geometry` | C11–C28; robust quantized PSD factorization and legal shape cliffs | revisit physical basis conditioning or cone dual atom-generation |
| **Delta: inverse radical arithmetic** | `goal-inv-radical` | I01–I08, I19–I22; fused-square, rational numerator, phase-aware middle products | register/phase lifetime superoptimization |

You may instead run independent generalist agents using this backlog, but require **different research portfolios/initial architectures** to avoid four agents converging on the same polished incumbent. Goal-score labels of 13/17/20 need not dictate which task an agent chooses; the geometric mean rewards whichever task supplies real gains.

### 7.4 Resource and workspace hygiene

- **Only write in your assigned experiment directory.** Do not change `from_host/`, accepted round-4/5 artifacts, shared Git branches/index, shared environment installations, or another active agent's files. Set experiment-local `TMPDIR`, caches, result locations, compilation intermediates. Reuse completed artifacts read-only, with attribution.
- Cap math-library threads and GPU jobs to coexist with other agents. The completed research uses `flock /tmp/r6-fullverify.lock` to serialize expensive complete verifications; coordinate via that lock (or confirm the actual existing agreed lock path), and use a small number of workers so overnight runs do not exhaust memory.
- The host counts **qubits physically present**; do not count live input words as free because their semantic labels changed, and don't subtract a work bank without an explicit reversible erasure schedule. DQI **Q is squared**, so scratch peaks are expensive.
- When creating new KMB/KMX oracles, independently test arbitrary incoming public operands, cleaned scratch, original classical mask/register preservation, phase and reverse execution. A function that works on zero-initialized outputs may still fail the actual XOR oracle contract.
- Use persistent artifact manifests including source provenance, license notices for adapted public code, fixed numerical seeds, full command-line arguments and compiler/driver/runtime versions.

### 7.5 Evidence ledger and honest stopping states

Maintain `RESULTS.md` throughout, not only at the end. For each attempt record:

```
Hypothesis and its mathematical invariant
Prior evidence / how this differs from previous failed designs
First decisive cheap test and result
Correctness proof or independent model validation
Measured native component T/Q, including forward and inverse
Projected complete T/Q/R or Lambda and why the assumptions are optimistic
Exact verifier PASS/REJECT, shot count, digest, salt
Original and new task cost, delta, contribution to combined score
New insight; why the result is local or general; reproducible next experiment
```

Use these terminal statuses:

- **Verified win:** frozen complete candidate passes the unchanged full native verifier, improves the task's score and passes independent validation. Provide an assembled archive only if the matching other components and packaging were checked.
- **Validated component breakthrough:** an algebraic invariant or reversible subroutine passes native tests, with actual resource counts, but no full task-cost win. Preserve it with a precise integration and break-even plan.
- **Scoped negative:** a **specified** architecture or exact dictionary is ruled out by measured costs, native failures or a genuinely applicable bound; state what different design would evade the failure.
- **Unresolved:** mathematical/numerical hypothesis remains underexplored; report missing test rather than embellishing feasibility or impossibility.

Do not select equivalent circuits or masks by favorable canonical-hash or random-salt outcomes. DQI errors affect `R`; use independent law-sampled failure estimation and independent full verification salts. Do not pool independent salts as a substitute for the host's per-run confidence factor.

### 7.6 Exploration policy when the scoreboard is flat

At a plateau, do one or more of the following rather than stopping:

1. **Change representation**, not a hyperparameter: DQI Euclid -> C-fraction/tangent/mask-functional; chemistry fixed tree -> new conditioned quantized projector dictionary or new QROM code; inverse Horner -> narrow fused square/rational phase.
2. **Transpose the computation**: push a classical output functional backward through a reversible polynomial map, so the phase is obtained without materializing the old intermediate output.
3. **Exploit exact task-specific promises**: reachable normalized inputs, squarefree/error weight distributions, one-body/two-body separate residuals, valid address care sets, or unused clean quantum bits—only after proving that the promise is respected coherently.
4. **Change where information lives**: read-only syndrome words can become cofactors; phase bank can borrow temporarily clean output bits; a measured bit can be erased with HMR only if its phase is correctly repaired.
5. **Cross-fertilize completed results**: e.g. care-set streamed erasure from `inv-dyadic` enabled `chem-qrom`; trim/retained-inverse logic may reduce `dqi-gap`; dual PSD certificates can guide **new atoms**, not merely reject old ones.
6. **Use a new algebraic identity**: small finite-field model and automated proof/reference first, full circuit only after a compatible liveness estimate.
7. **Try a representation that deliberately sacrifices a little numerical accuracy or Toffoli count to save qubits**, or vice versa. Track the actual product and its discontinuities, not one metric.
8. **Actively search for counterexamples to your own pessimistic resource bound.** Identify hidden assumptions: fixed rank/C, fixed bit precision, fixed tree, same error rate, materialized intermediate, perfect decoding, fixed oracle batch, and no code reuse.

### 7.7 Final integration and handoff

At the end of an unattended run, leave **one short status paragraph** at the top of `RESULTS.md` listing best verified candidate, task cost, combined score if integrated, and exact paths. Preserve the selected incumbent as an immutable fallback. If a combined archive is produced, independently reopen it, check member/digest integrity, run the unchanged complete metric including all default shots, and report its **actual** combined score. Do not claim a score of 13/17/20 from multiplying component projections.

---

## 8. Most useful references and provenance

### Appendix A — completed local research to read first

Source reports supplied by the user on October 24, 2026; repository names below are the paths used **within** those reports. If a path has since changed, find it from the repository source manifests and hashes rather than inferring it from this Markdown document.

**DQI**

- `initial_experiments/5/dqi-gap/RESULTS.md`: currently best verified complete DQI; schedule, gap-aware arithmetic, 1004Q liveness, zero-error native acceptance.
- `initial_experiments/5/dqi-repair/RESULTS.md`: exhaustive diagnosis of degree-two quotients and first genuine reversible exception block.
- `initial_experiments/5/dqi-core/RESULTS.md` and `initial_experiments/5/dqi-core/DESIGN_MEMO.md`: accepted 961Q trimmed approximate decoder, 25,596T saved, exact middle-product attribution and failing ratio designs.
- `initial_experiments/5/dqi-core-phase/RESULTS.md`, `DERIVATIONS.md`: Riccati, C-fraction, Wronskian H, 46T point kernel, component and census data.
- `initial_experiments/5/dqi-census/RESULTS.md`, `PHASE1.md`, `BUILD_PLAN.md`: detailed candidate cost budgets, fixed-degree modular inversion base cases, algebraic identities and outstanding primitives.
- `initial_experiments/5/dqi-phase/RESULTS.md`, `dqi-hybrid/RESULTS.md`: negative workspace results for separate Psi/F and literal BM banks; exact liveness controls.
- `initial_experiments/5/dqi-dispatch/RESULTS.md`, `RESULTS.matrix-history.md` if present: classical general-quotient transfer and losing separate-chart/MAC implementations; exact field/head component circuits.
- Round 3 DQI stream decoder, round 4 Codex compact decoder, and pinned `initial_experiments/3/vendor/uncompute` under Apache 2.0.

**Chemistry**

- `initial_experiments/5/chem-census/RESULTS.md`: currently best reported structured-banked full chemistry PASS at 696.943M; rank/width jumps and all invalid continuation/growth controls.
- `initial_experiments/5/chem-qrom/RESULTS.md`: 4-row nonlinear decoder, bank-bit elimination, prefix HMR, row/permutation experiments and two-salt native compiler win.
- `initial_experiments/5/chem-audit/RESULTS.md`: fixed-shift **universal** normalization lower bound, banked incumbent margin refinement, and assembled archive with *older trimmed DQI*.
- `initial_experiments/5/chem-discrete/RESULTS.md`, `chem-psd/RESULTS.md`: exact-word/spectrum primal-dual PSD certificates, refitting, numerical basins, rank/column frontiers and low-Lambda failures.
- `initial_experiments/5/chem-tree/RESULTS.md`: legal rooted tree semantics, 17 binary lattices, re-encoding, conditioning/singular-direction analysis and nine-bit width cliff.
- `initial_experiments/5/chem-qfit/RESULTS.md`, `chem-blocks/RESULTS.md`, `chem-blocks/jad_search/RESULTS.md`: accepted-basin calibration, continuous-vs-discrete fitting and resource frontiers.
- Round 4 Codex selected `norm_cycle2_signed`, physics/gauge/oracle builders and `NORMALIZATION_BOUND.md`.

**Inverse square root**

- `initial_experiments/5/inv-census/RESULTS.md`, `BUILD_PLAN.md`: exhaustive passing **narrow fused-square** *numeric* formula, injection-based exponent marker, no winning complete circuit yet.
- `initial_experiments/5/inv-phase/RESULTS.md`: exact-rational quadratic with much larger synthesis margin; verified but losing signed-row/split phasing.
- `initial_experiments/5/inv-dyadic/RESULTS.md`: independent reproduction of existing prefix table, signed overflow controls and streaming ANF cleanup.
- `initial_experiments/5/inv-width/RESULTS.md`: exact 90Q liveness at first final Horner row, failed low-ancilla/segmented designs.
- Round 4 Codex `zero_head_affine.kmb` and its builder, original exhaustive integer phase reference.

### Appendix B — public technical literature / implementations

These are external context, not measured competition scores. Read a paper's specific subroutine and source code rather than adopting its whole-task cost estimate.

1. Gu, Kemeny, Yelin and Lukin, **Discovering state-of-the-art quantum circuits with large language models** (2026), https://arxiv.org/abs/2610.22367; full examples and Apache-2.0 sources at https://github.com/andigu/uncompute . Includes RS DQI, Fe2S2 lookup compression, HHL-angle/coulomb approximation and branch-history reconstruction.
2. Khattar et al., **Leveraging verifiable structure in quantum circuits to find more efficient compilations with agentic AI** (2026), https://arxiv.org/abs/2610.23832 . Host task framing, kickmix semantics and expected verification behavior.
3. Khattar et al., **Verifiable Quantum Advantage via Optimized DQI Circuits** (2025), https://arxiv.org/abs/2510.10967 . Implicit Bezout coefficients and synchronized Euclidean designs; understand what is and is not a theoretical space bound.
4. Low et al., **Faster electronic structure quantum simulation by spectrum amplification** (2025), https://arxiv.org/abs/2502.15882 . DFTHC, SOS normalization and cost-model context; published optimizable shifts are not necessarily allowed in this fixed-shift competition.
5. Reed–Solomon decoder algorithms and partial approximants, https://fredrikj.net/fdoc/html/tutorials/algebra/reed-solomon.html . Reference classical algorithm families, not necessarily reversible winners.
6. AlphaTensor-Quantum field-multiplication circuit optimization, https://www.nature.com/articles/s42256-025-01001-1 . Inspiration for specialized bilinear circuit search, not a verified 12-bit GF multiplier for this task.
7. Qualtran/kickmix and exact verifiable circuit model, https://github.com/quantumlib/Qualtran . Use the bundled pinned host implementation as authoritative for scoring.
8. Published LLM-research-process discussion, https://www.anthropic.com/research/designing-efficient-quantum-circuits . Illustrates why **representation and uncomputation**, rather than only gates, can deliver large long-horizon gains; the paper/repo provide reproducible methods.
9. Google Quantum AI competition overview and discussion, https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation . New public notebooks/discussions may reveal different algorithms; validate sources and license terms independently.

---

## 9. Optional goal-mode instruction to paste before referencing this backlog

> You are conducting a **long-horizon, autonomous research campaign** on the Verifiable Quantum Circuit Compilation competition. Read `quantum_goal_mode_research_backlog_2026-10-24.md` and the completed experiment reports it references. Aim for local combined score **TARGET** (13, 17, or 20), but make honest, score-backed progress even if you cannot hit the target. Begin from the best *independently verified* per-task components and assemble a full metric-verified baseline. Explore **multiple distinct mathematical architectures**, prioritizing unbuilt leads with clear resource tests. Maintain a checkpointed research ledger and a frozen best submission. Do not declare the goal impossible based on failure of a single representation, optimizer or cost projection; explicitly challenge those assumptions and pivot to another idea. Equally, do not claim a mathematical, resource or score win without the appropriate proofs/tests and full native verification. Work only in your assigned directory; read other **completed** work read-only and do not inspect active agents. Keep CPU/GPU usage compatible with other agents; serialize full verifier runs under the agreed lock. Preserve important intermediate components and candid negative evidence. Provide a concise final status and actual verified best task costs, combined score, source hashes and reproduction commands. Continue productive research independently rather than stopping after the first negative experiment.

**Final reminder:** The leader's reported ~18.9, or a hypothetical future 20, disproves *nothing specific* about a proposed mechanism; conversely, all current local failures are insufficient evidence for a general impossibility theorem. High-confidence progress comes from **new reversible representations, honest complete objective accounting, diverse searches, and frozen verified artifacts**.