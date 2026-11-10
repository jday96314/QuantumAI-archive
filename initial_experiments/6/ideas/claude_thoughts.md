# Round-6 /goal playbook: approaches worth exploring

Prepared 2026-10-24 from all round-5/6 reports (Claude- and ChatGPT-directed agents). Intended for autonomous agents. Read sections 1–3 before starting.

## 1. Where we stand

Best verified pieces (all pass full-default native verification):

| Task | Best verified artifact | Cost | Gain vs round 4 | Resources |
|---|---|---:|---:|---|
| Inverse | round-4 Codex `zero_head_affine.kmb` (SHA `27b14e54…`) | 182,956.53 | 1.00× | Q 90, mean T 1,782, 34 rotations |
| DQI | gap-aware decoder `gap_best_n70_m4095.kmb` (SHA `8a0675d5…`, built by `build_best.py`) | 1.4359e12 | 1.113× | Q 1004, mean T 1,195,802, **0 failures** in three 102,400-shot runs, R 1.1913 |
| Chemistry | chem-census `candidates/structured` (digest `62f00bbf…`) | 696,942,847 | 1.017× | R1/B719/C350 banked words, J4 structured ANGLE compiler; T 31,608, Q 454, λ_eff 48.567 |

These three combine to **≈9.76** (round-4 reference 9.3677). `chem-audit/submission.zip` (9.47 locally) uses the older trimmed DQI decoder (Q 961, 1.560e12) and the banked chemistry (708.16M). It is therefore not the best available combination; reassemble it first (X1).

Find the gap-aware decoder with `find initial_experiments/5 -name gap_best_n70_m4095.kmb`.

The gap-aware decoder decodes 2–3 per million shots wrongly, about 0.25 expected per 102,400-shot run. Its R normally sits at the zero-failure floor (1.1913), so DQI cost is now close to deterministic and small DQI gains are measurable without salt noise. A rare failure on the private salt raises R slightly.

## 2. What 13, 17 and 20 require

Required product of task improvements over the current best pieces: **13 → 2.36×, 17 → 5.28×, 20 → 8.60×.**

One balanced split for each target:

| Overall | DQI cost (factor) | Chemistry cost (factor) | Inverse cost (factor) |
|---|---:|---:|---:|
| 13 | ≤0.72e12 (2.0×) | ≤657M (1.06×) | ≤163k (1.12×) |
| 17 | ≤0.55e12 (2.6×) | ≤498M (1.40×) | ≤126k (1.45×) |
| 20 | ≤0.48e12 (3.0×) | ≤422M (1.65×) | ≤105k (1.75×) |

What those numbers mean physically:

| Target | Example operating point |
|---|---|
| DQI ≤0.72e12 | ~745k mean T at 900 Q, or ~780k at 880 Q, with zero failures |
| Chemistry ≤498M | T×Q ≈10.3M at λ_eff 48.6, i.e. a B≤511-class table |
| Inverse ≤126k | e.g. Q 70 × 1,800 Toffoli-equivalents |

**Honest read.** 13 needs one architecture-level DQI win (~2×) plus modest gains elsewhere; sections 4–6 contain two credible routes to that. 17 and 20 need breakthroughs in all three tasks. The approaches marked "high upside" are the candidates, and none has a demonstrated path yet. The public leaderboard (18.9) shows such gains exist somewhere.

**Goal-setting:** give each single-task agent a task-cost goal from this table, not an overall score.

## 3. Rules for every agent

- **Files.** Write only inside your assigned directory. Read anything else read-only. Never modify verifiers, confidence or scoring code, or any `from_host/` file.
- **Objectives.** Score with the actual formulas:
  - DQI: `R(eps_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson).
  - Chemistry: `lambda_eff * T_total * Q_total`.
  - Inverse: `Q * (mean_Toffolis + total_T/4)`, with rotation synthesis priced from the exhaustive minimum margin.
- **Acceptance.** Accept only on full-default native verification plus one independent predeclared salt:
  - DQI: 102,400 shots.
  - Inverse: 200,540,928 shots.
  - Chemistry: native ANGLE recovery and physical checks.
  - Never search salts or hashes, and never select a favorable salt.
- **Shared machine.** 8 cores and one GPU. Respect your CPU cap and wrap full verifier runs in `flock /tmp/r6-fullverify.lock`. Use the GPU one job at a time.
- **Charge everything.** Every emitted branch counts (no frequency discount for rare cases). Report width at the true peak, with a liveness table. Use exhaustive synthesis margins, not sampled ones.
- **Persistence.** Do not stop because a first draft loses to a multi-round incumbent. Drop a line only when its key uncertainty is resolved against it, and record that negative result scoped to the design tested.
- **Reporting.**
  - Report a numeric first checkpoint early, and keep `RESULTS.md` current.
  - When a verified improvement lands, put `ACCEPTED: <artifact path> <SHA>` at the top of `RESULTS.md`.
  - Attribute reused Apache-2.0 code (andigu/uncompute rev `64b5a13a…`; host verifiers).

## 4. DQI approaches

**Exchange rates near the current point:**
- 10 qubits ≈ 2% of cost.
- 12k Toffolis ≈ 1%.
- The syndrome register alone is 840 Q.

| ID | Approach | Potential on DQI | Risk | First milestone |
|---|---|---:|---|---|
| D1 | Polish the gap-aware decoder | 1.10–1.25× | Low–medium | Gap-aware rebuilt with core trims, verified |
| D2 | **Blocked C-fraction decoder** (high upside) | 1.5–2.2× | Medium–high | Classical blocked model with 0 host failures + one-block profile |
| D3 | **Resident σ/ω in syndrome storage** (high upside) | up to ~2.4× | High | Alias-safe update proven with ≤~100 extra Q |
| D4 | Scan-body reductions | 1.03–1.15× | Low–medium | AND-count search result for the body |
| D5 | Width-first trims on any core | 1.05–1.10× | Low | Liveness table at the peak |
| D6 | Arithmetic kernels | 1.02–1.08× | Medium | Better middle product at a used size |

### D1. Polish the gap-aware decoder

Current: Q 1004, mean T 1,195,802 (core 824,740 + scan 371,150). Schedule: block-3 at k=1,4,7,10; two-slot at 13,15,17; one-pass singles at 19–34.

The peak (k=10 Schur head) is 1004 qubits:

| Register | Qubits |
|---|---:|
| Syndrome | 840 |
| Copied head | 84 |
| Field windows | 36 |
| History | 26 |
| Carry | 12 |
| Emitter | 3 |
| Controls | 2 |
| Chart control | 1 |

Its k=10 block costs 43,470 T, versus 34,953 for the old regular block. Cofactor routing costs 14,049 versus 8,982 T.

Moves:

1. **Port dqi-core's trims.** Trimmed middle products and budgeted inverse retention saved 25,596 T on the old core (`dqi-core/components/trimmed_middle.py`, `budget_inverse.py`).
2. **Remove or shrink the 84-bit copied head.** One-pass reconstruction already avoids it in late slots. Rerun the schedule dynamic program with head-free or shorter-head early blocks. Going from 1004 to ~965 Q is worth ~8%.
3. **Compress the 26 history bits.** dqi-dispatch's fixed-schedule contract uses 21 metadata bits (6 count + 2×6 gap tags + 3 mode) and had 0 wrong decodes in 2.1M classical cases. Its completion-counter primitive saved 13 Q at the old peak.
4. **Cheapen routing with dispatch's measured components:** shared cofactor with degree-≤3 multipliers (3,858 T at d19), routed Schur head (2,454 T), carry-finishing slot (1,572 T) and merged prefix realignment (684 T). Also use gap-aware's own Laurent cofactor.
5. **Add the 70-point interpolation skip** in the scan (~−6k T, zero qubits).

Goal: ≤1.25e12 (1.15× over current); stretch ≤1.15e12.

Pitfall: rerun the carried census and the residual-path checks after every change, so zero host failures are preserved.

### D2. Blocked C-fraction decoder (high upside, evidence-backed)

The phase-aware study built a reversible C-fraction decoder:
- Its core lives in **878 Q including the 840 syndrome qubits**.
- One coefficient word per peeled moment is stored in the syndrome array.
- B and Ω are synthesized in place (57,664 T per direction), followed by an 80 T/point Forney scan.
- Native small oracles (N6, N14) pass on supported charts.

The problem is cost: serial extraction is 503,655 T per direction, giving ~1.45M T at ~884 Q.

What it needs:
- **Blocked extraction.** Peel several coefficients per block from a few head words (as the compact Euclid's six-moment heads do), then update the tail once per block instead of once per peel.
  - Extraction ≤292,699 T per direction gives **1.5×** over the current best.
  - Extraction ≤164,158 gives **2×** (both at ~884 Q with zero failures).
- **Chart handling for zero heads and degree gaps.** Generic charts currently fail on ~1.3% of shots (253–291 wrong pairs in 20,000). Reuse dqi-dispatch's fixed-schedule local-transfer algebra (handles (2,1), (1,2), (3), deferrals and completion; 0 failures in 2.1M), or gap-aware's carry charts.
- **Width:** stay at or below ~900 Q.

First milestone: a classical model of blocked extraction with charts showing 0 failures on the 102,400 host cases (≤~2 per million on law streams), plus an exact T/width profile of one block at k≈10.

Start from the phase-aware report (`cfraction_decoder.py`, `DERIVATIONS.md`), dispatch's `local_transfer.cpp`, and the round-4 compact-Euclid block heads.

Risks: a per-block tail update may cost as much as a full series division, and chart handling may add width.

### D3. Resident σ/ω in syndrome storage (dqi-census "C", high upside)

- **Target:** the census conditional target is Q 940, T ~572k, giving 0.602e12 (**2.4×** over current).
- **Missing primitive:** an alias-safe update in which σ and ω live in the processed syndrome words with no separate τ bank. At a generic 2k boundary, k σ words and k ω words fill the 2k processed words, so a full τ has nowhere to live.
- **Rejected attempt:** the straightforward reference realization needs 1,707 Q.
- **Possibly useful:** the Riccati tangent identity (DB = Ω, DΩ = 0, from the phase-aware `DERIVATIONS.md`) might let τ be reconstructed instead of stored.

First milestone: a classical proof plus a small native demonstration of the alias-safe update using ≤~100 extra qubits. If the primitive exists, D3 can share chart handling with D2.

### D4. Scan-body reductions (apply to whichever core wins)

- **Body AND-count search.** Today's bodies cost 88 T (inverse-product) and 80 T (Forney). Gu et al.'s Hermitian GF(64) scan reaches 33 T/point by inverting in the subfield F_8. Search GF(4096) ⊃ GF(64) ⊃ GF(8) towers for the gated inverse-trace, which costs 46 T today. Each −8 T/point ≈ −33k T.
- **46 T inverse-only body with Wronskian H.** Built and verified at 46 T/66 Q. Preparing H costs 100,080 T (build + erase) at 897 Q. It pays only if H replaces a polynomial that is already stored, not as an extra bank.
- **11 T F′ scan.** Needs an in-place degree-35 map A → (A·B′)⁻¹ mod B with ≤~74 extra qubits; the degree-2 base case exists (960 T/69 Q). Break-even is ~155–158k T per direction at unchanged width. Fold β^70 into the classical masks.
- **70-point interpolation skip** (classical masks, zero qubits).

### D5. Width-first trims (any core)

Profile liveness at the true peak first. Then:
- Recompute completion state instead of storing flags.
- Reuse cleared flags as scan scratch (gap-aware already does this).
- Keep scale words in syndrome words known to be zero at that stage.

### D6. Arithmetic kernels

Series and middle products are ~50% of compact cores; a field multiply costs 42 T and an inverse 66 T. Search better transposed or middle-product algorithms at the block sizes actually used. Expect small percentage gains everywhere rather than one large win.

### DQI: closed or rejected (don't retry without a new idea)

- **Explicit-matrix/Bennett-bridge dispatch:** ≥1,107 Q.
- **Separate moment charts per carry case:** +326k T.
- **Full running-ratio recurrence:** 3–11× worse per block.
- **Block-4 in the old frame:** 984 Q (may become feasible if D1/D5 free width).
- **Mechanical reverse asymmetry:** 0 T saving.
- **Ψ/F transplant as built:** 1,383–2,931 Q.
- **Frobenius F′:** 11,751 Q.
- **Whole-bank F′ XOR/swap at degree 35:** +420 Q.
- **Wide BM/Ψ:** 1,837 Q.
- **Lazy-BM handoff with zero-cost splice:** ≤1.23×.
- **General Hankel/PGZ:** 2,200 Q.
- **Streaming residue:** ≥12.7M T.
- **Retaining a separate H bank:** loses on Q².

## 5. Chemistry approaches

**Fixed facts:**
- λ_eff ≥ 46.786 is proven for every factorization at the fixed shift, so normalization has ≤3.8% left.
- R·C ≥ 290 is necessary at 20 mHa.
- Inner preparation ≈ R·C·⌈(B+1)/16⌉ is ~55% of T and moves only with shape.
- One-body and two-body error budgets are independent.

**Exchange rates:**
- 300 T ≈ 1%.
- 4.5 Q ≈ 1%.

| ID | Approach | Potential | Risk | First milestone |
|---|---|---:|---|---|
| C1 | Rank trim C350→C340 at fixed words | 1.01–1.02× | Low | Native pass at C340 |
| C2 | ANGLE compiler improvements | 1.03–1.10× | Medium | ≥600 T saved at J4 |
| C3 | **Λ-aware deletion continuation to B≤511** (high upside) | 1.15–1.35× | Medium–high | Passing state at a 16-row threshold with λ_eff ≤56 |
| C4 | **Lattice-native well-conditioned dictionaries** (high upside) | up to ~1.4× | High | Exact-word dictionary with PSD dual bound <20 mHa at T×Q ≤12M |
| C5 | Run the prepared robustness branch | Unknown | Medium | First robust refit from an incumbent seed |
| C6 | Steer by the host-threshold map | Enabler | Low | — |

### C1. Rank trim at fixed words

chem-census measured the PSD rank-C capacity of the unchanged B719 words:

| C | Capacity error (mHa) | Passes 20 mHa? |
|---:|---:|---|
| 340 | 19.910 | Yes |
| 320 | 21.36 | No |

Rebuild at C340 (or the smallest passing C after integer word moves) with the structured compiler, then re-optimize weights and offsets for λ_eff under the error cap. Inner preparation falls 45 T per column, ~450 T (1.4%) at C340, provided λ_eff holds. This takes hours.

### C2. ANGLE compiler

- **Current:** chem-qrom's J4 hybrid costs 309 T per oracle (292 for oracle 0 with feature-aligned rows) across 19 oracles, 11,708 T over both passes. A 1.10× gain needs ~2.4k fewer T.
- **Joint optimization:** choose quadratic decoder factors, omitted columns and a shared address order jointly across all oracles.
- **Per-oracle bases:** extend oracle 0's feature-aligned trick to every oracle.
- **Bitwise or linear factoring:** chem-audit's zero-saving certificate covers only whole-payload equality, so these remain open.
- **Re-tune J after any dictionary change:** persistent word storage is a width lever (J4 saved 24 Q).

### C3. Λ-aware deletion continuation toward B≤511 (top chemistry lever)

Start from the incumbent integer basin and repeat:
1. Delete rows in small batches.
2. Repair words.
3. Refit rank-C weights.
4. Periodically run Λ-constrained weight and offset refits.

Evidence so far:
- At B703, 240 repair rounds reached 19.98 mHa; the constrained refit there passed physics at λ_eff 56.4.
- Sequential deletion without repair degrades quickly (B639: 56.6 mHa).

Where the payoff sits:
- Each 16-row threshold cuts ~560 T (351 T of inner preparation plus ANGLE).
- B513→512→511 drops Q by ~18 (478→460 with J7).
- At B511 (T ≈23.4k, Q ≈460 with J7, lower with J4), break-even λ_eff against the current best is ≈64.6. λ_eff ≤56 gives ≈1.15×; λ_eff ≈48.6 gives ≈1.33×.

Intermediate states rarely beat the incumbent, so aim at B≤511 directly and judge progress by the trajectory of (error, λ_eff). Allow re-adding residual-directed atoms and incumbent atoms as repair moves. R2/B351/C175 (conditional 1.45× at incumbent λ) is a second endpoint worth tracking.

### C4. Lattice-native, well-conditioned dictionaries

R2/B351/C189 passes in continuous form (19.98 mHa, conditional 1.40×), but rounding destroys it:
- Its feature blocks are badly conditioned (condition numbers 1,525 and 2,210, versus 172 for the incumbent).
- The rounded span bound is 140 mHa; the best integer repair reached 54.3 mHa.

Build dictionaries directly on the 8-bit lattice:
- Penalize small singular directions and Λ.
- Use tensor-Huffman trees, which minimize isotropic quantization sensitivity.
- Use block-PSD dual certificates to discard dead dictionaries early.

### C5. Run the prepared robustness branch

chem-audit's `run_robust.sh`, `integer_refine.py` and `integer_grow.py` were never run. Run them seeded from the incumbent or from C3 states, not from scratch.

### C6. Host-threshold map

`chem-census/results/thresholds.json` records every B479–719 point and change. Avoid 9-bit angles (the k2 cliff takes temporary width from 321 to 640) and more than 257 vectors per block.

### Chemistry: closed or rejected

- **Scratch float fits:** 28–47 mHa before rounding.
- **Unregularized multiblock fits:** Λ ≈ 1,000+.
- **Exact prior dictionaries:** B180, B224, B351, B479 and B479-repair are excluded by PSD or span certificates.
- **Growth from small integer dictionaries:** excluded by span bounds at 240 and 35 mHa.
- **Subtree pruning or row permutations:** 0 T saved.
- **Hamming or lexicographic orders, 128 Givens gauges:** no gain.
- **Direct ANF compiler:** loses.
- **Uniform one-body rounding:** 123 mHa.
- **Nine-bit shifted grid:** k2 width cliff.
- **C≤320 with unchanged words:** fails accuracy.
- **Any R·C<290:** fails the necessary rank condition.

## 6. Inverse-square-root approaches

**Facts.** The 90-qubit peak is in the first radix-4 row of the final Horner product:

| Live register | Qubits |
|---|---:|
| Input | 22 |
| Constant remainder | 23 |
| Slope | 17 |
| Output headroom | 10 |
| Exponent flags | 5 |
| Curvature | 4 |
| Reclaimed | 3 |
| Carries | 3 |
| Counters | 2 |
| Booth temporary | 1 |

| Stage | Cost (Toffoli-eq) |
|---|---:|
| Normalization | 210 |
| Lookup + erase | 304 |
| Horner | 1,046 |
| Output shift | 217 |
| Affine restore | 16 |
| 34 rotations | 251 |

**Exchange rates:**
- 1 qubit ≈ 1.1%.
- 20 T ≈ 1%.

| ID | Approach | Potential | Risk |
|---|---|---:|---|
| I1 | **Fused narrow quadratic correction** (high upside) | 1.1–1.3× | Medium |
| I2 | In-place exponent encoding (frees 5 flag qubits) | ~1.05× | Medium |
| I3 | Spend the exact-rational margin on narrower registers | 1.02–1.08× | Medium |
| I4 | Rerun the lifetime search after I1–I3 | Enabler | Low |
| I5 | Coarse-coordinate higher-order corrections | 1.05–1.15× | Medium |
| I6 | Normalization and mean-T accounting | 1.02–1.05× | Low |

### I1. Fused narrow quadratic (inv-census, numerically verified, not built)

```
u   = floor((z + 16) / 32)        # 9-bit unsigned
h   = round(u^2 / 64)             # 11-bit
num = c0*2^23 - c1*z*2^11 + c2*h*2^16
phase_units = num * 2^(s+4) mod 2^64
```

- This passes on all 3,133,452 inputs with error/tolerance 0.805 and a margin of 137.9M units; a q=4 variant reaches 0.696 and 237.2M.
- It replaces the curvature Horner stage with a ~9-bit square and a narrow product.
- The slope then comes straight from the table, so the 17-bit updated slope need not stay live across a second product.
- Order the work to minimize the peak: accumulate c2·h first, erase c2 and h, then accumulate c1·z.
- inv-census's `BUILD_PLAN.md` lists the kernels and the coefficient-liveness problem.
- Budget: Q 80 with T 1,400 gives ~140k (1.30×); a realistic Q 84 with ~1,750 T gives ~1.07–1.1×.

### I2. In-place exponent encoding

- **Map:** E = M for k=0, else E = M − 2^21 + 2^(k−1). It is injective over the domain, and k = ctz(E)+1.
- **Saving:** it frees the 5 flag qubits that are live at the peak.
- **What to price:** the on-demand decode, and the change it causes to the local coordinate.

### I3. Spend the exact-rational margin

Evaluating the quadratic exactly (no per-row rounding) gives 14.77× more synthesis margin (295.3M units). Spend it on:
- narrower slope and accumulator registers at the peak,
- fewer output bits and rotations,
- or fewer table leaves.

Do not just bank it as synthesis savings; that is worth only ~34 units.

### I4. Rerun the lifetime search

After I1–I3 change the arithmetic, rerun the schedule search with these known tools:
- streamed ANF erasure (31 instead of 83 scratch qubits at an 11-bit address, same Toffolis),
- the circular, wrap-corrected result (no output headroom; 278 T kernel, 37 rotations),
- delayed constant bits,
- conditional carry repair.

Headroom removal alone failed before because arithmetic scratch set the peak; with new arithmetic, recheck.

### I5. Coarse-coordinate higher-order corrections

Small corrections evaluated on a coarse coordinate are cheap, as I1 shows. Use them to cut table leaves (the lookup costs 267 T; the incumbent table has 269 binary-prefix leaves) and coefficient widths together.

### I6. Normalization and mean-T accounting

- Normalization costs 210 T over both directions; the reverse currently runs at 101 static versus 90 mean.
- Look for measurement-based uncomputation where fix-ups execute only on half of the outcomes.
- Audit every uncompute for this pattern.

### Inverse: closed or rejected

- **Incumbent schedule grids:** 423 + 362 + 189 schedules; best 88 Q at 183.2k.
- **Segmented phasing:** the carry spans 6–9 bits.
- **Dirty-input borrowing.**
- **Bit-serial Horner:** 448k.
- **Table-linear:** 2k+ T just to load.
- **Direct quadratic with separate term lifetimes:** 2.0M.
- **Walsh exponent rotations:** 307 rotations.
- **Split-bank phasing:** 84 Q, 231k.
- **Direct row phasing:** 80 Q, 318k.
- **Borrowed or circular headroom alone:** 190–196k.
- **Adaptive dyadic local coordinates:** 92 Q.
- **Mixed degree.**
- **26-bit outputs:** fail accuracy.

## 7. Cross-cutting

- **X1. Assemble the best archive now, and again after every acceptance.**
  - Combine gap-aware DQI, chem-census `structured` chemistry and the incumbent inverse, then run `from_host/metric.py` (≈9.76 expected).
  - Tooling: round-4 Codex `package_submission.py` and `audit_submission.py`, plus chem-audit's `assemble_archive.py` and `package_metric.py`.
- **X2. Recheck failure models after any DQI change**, using dispatch's carried census, gap-aware's `carried_census.cpp` and its residual extraction.

## 8. Suggested overnight split (5 agents)

| Agent | Directory | Work | Goal | CPU / GPU |
|---|---|---|---|---|
| 1 | `dqi-polish` | D1 + D4 small items + D5 | DQI ≤1.25e12 verified | 2 |
| 2 | `dqi-cfrac` | D2, falling back to D3 | DQI ≤0.72e12 | 2 |
| 3 | `chem-continue` | C1 first, then C3 (C4/C5 as alternatives) | Chem ≤657M; stretch ≤600M | 1 + GPU |
| 4 | `chem-compile` | C2 | Chem ANGLE −2.4k T | 1 |
| 5 | `inv-fused` | I1–I4 | Inverse ≤163k; stretch ≤140k | 2 |

- Whoever finishes first runs X1.
- Stagger full verifications under the lock.
- If you have capacity for a sixth agent, give D3 its own directory.

## 9. Report index

Find reports with `grep -l "<title>" initial_experiments/5/*/RESULTS.md` where the directory is unknown.

| Area | Directory or report title | Key content |
|---|---|---|
| DQI | "Round 5 gap-aware compact DQI decoder" | Current best decoder, schedule DP, residual-path checks |
| DQI | "Phase-aware decoding research" | C-fraction decoder, Wronskian H, 46 T body, DERIVATIONS |
| DQI | `dqi-census` | Eight priced architectures, resident σ/ω target, F′ degree-2 map |
| DQI | `dqi-dispatch` | Fixed-schedule classical model (0 failures), shared components, completion counter |
| DQI | `dqi-core` | Trims, attribution, accepted trimmed decoder |
| DQI | `dqi-repair`, `dqi-phase`, `dqi-hybrid` | Earlier diagnosis, Ψ/F costs, liveness |
| Chemistry | `chem-census` | Current best chemistry, threshold map, B703 continuation, rank audit |
| Chemistry | `chem-qrom` | Structured four-row quadratic ANGLE compiler |
| Chemistry | `chem-audit` | Universal λ bound, banked fit, subtree certificate, packaging scripts |
| Chemistry | `chem-discrete` | B711 controls, PSD certificates for B479 |
| Chemistry | `chem-psd` | R2/B351/C189 continuous pass, resource table |
| Chemistry | "Shared-tree chemistry representation research" | Tree topology, tensor-Huffman trees, conditioning diagnostics |
| Chemistry | `chem-qfit`, `chem-blocks` | Earlier fits and host-exact frontiers |
| Inverse | `inv-census` | Fused-square model, in-place exponent map, BUILD_PLAN |
| Inverse | `inv-dyadic` | Dyadic fits, streamed ANF erasure |
| Inverse | "Phase-oriented inverse square root" | Exact rational quadratic margin, row-phase identities |
| Inverse | `inv-width` | Incumbent liveness, schedule grids |

External references:
- Gu et al., arXiv:2610.22367; repository https://github.com/andigu/uncompute.
- Khattar et al., arXiv:2510.10967.
- The competition paper, arXiv:2610.23832.