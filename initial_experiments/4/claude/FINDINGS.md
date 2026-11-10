# Experiment 4 (claude): what was tried and what happened

## Outcome

No better submission was produced. `submission.zip` is still experiment 3's verified
components, with an overall score of **9.0096**. Two goals were set during the round:
- 14.0, later cleared;
- 9.4, never reached.

The same round's Codex workspace did reach a verified **9.3677** (see `../REPORT.md`).
That disproved several limits I stated to the user and wrote into these docs. Each wrong claim
is recorded in **[MISTAKES.md](MISTAKES.md)**, quoting my original wording. The
[Corrections](#corrections-in-light-of-codexs-results) section below summarizes them.

Each item below is labelled with how strong its evidence is:
- **Measured**: a run in this workspace produced the number; the log or data file is named.
- **Host-verified**: the competition verifier produced the number.
- **Estimate**: arithmetic or a cost-model calculation only; no circuit or fit was built.

## Starting point

Experiment 3's verified costs, with each ratio taken against the competition baseline:

| Task | Score formula | Cost | Improvement |
| --- | --- | ---: | ---: |
| invsqrt | Q × (Toffoli + T/4) | 94 × 2,028.85 = 190,712 | 10.89× |
| DQI | R · Q² · T | 1.223 × 1,429² × 645,755 = 1.613e12 | 15.15× |
| Chemistry | λ_eff · T · Q | 48.895 × 29,440 × 526 = 7.572e8 | 4.43× |

The overall score is the geometric mean of the three ratios. Reaching 14.0 needs their product
to rise 3.75×; reaching 9.4 needs 1.136×.

## Chemistry

Almost all of the compute in this round went into chemistry.

### Cost model and sanity checks

- `costscan.py` reproduces the verifier's heuristic state-preparation cost and the
  angle-circuit model 2·K·(rows − 1), with K = ⌈75/J⌉.
- **Host-verified:** rebuilding experiment 3's chemistry submission with the copied pipeline
  (`chem/build_global_chem.py` and its helpers) reproduced 757,163,489 exactly.
- **Host-verified (resources only):** an end-to-end `chem/pipeline.sh` run on an unconverged
  R=2, B=351, C=175 plumbing fit (`q_plumb`, error 0.24) gave T = 20,793 and Q = 508.
  Experiment 3's shape (R=1, B=719, C=350) has T = 29,440 and Q = 526. That build failed
  accuracy as expected; it only checked the pipeline and the resource accounting, and its
  output directory was not kept.

### Two-block shapes: the main effort

**The idea.** The inner state preparation costs R·C·⌈(B+1)/16⌉. Two blocks with B=351 and
C=175 hold about as many parameters as R=1 with B=719 and C=350, but cost about 1.47× less
(estimate, consistent with the resource check above). The 20 mHa error cap decides whether
it is usable.

**Float fits from scratch: measured, all at or above about 0.021.**

| Run (log) | Method | Shape (R, B, C) | Float error₂ |
| --- | --- | --- | ---: |
| `L_r1b719c350` | L-BFGS with exact w projection (`lbfit.py`) | 1, 719, 350 | 0.0210 |
| `L_r2b351c175` | same | 2, 351, 175 | 0.0211 |
| `L_r4b180c88` | same | 4, 180, 88 | 0.0470 |
| `c_r2b360c175` | eigenvector init, then Adam (`init_eig.py`, `shapefit.py`) | 2, 360, 175 | 0.0280 |
| `c_r2b319c175` | same | 2, 319, 175 | 0.0391 |
| `c_r3b239c117` | same | 3, 239, 117 | 0.0360 |
| `c_r4b180c128` | same | 4, 180, 128 | 0.0375 |
| `bvpA` | variable-projection Adam, rank annealed 351→175 (`blockvp.py`) | 2, 351, 175 | 0.0215 |
| `bvpC190` | same | 2, 351, 190 | 0.0203 |
| `F_r2b351` | same, full rank | 2, 351, 351 | 0.0172 |
| `vp_r2` | variable-projection L-BFGS (`vpfit.py`) | 2, 351, 175 | 0.0305, stalled |
| `s_r1b719c350`, `r1b719c350` | Adam from random init | 1, 719, 350 | 0.072–0.075 |

`isdf_r2.npz` (greedy ISDF dictionary init, `init_isdf.py`) did no better than the
eigenvector init. The conclusions:
- At equal effort, R=2 matches R=1 (0.0211 vs 0.0210).
- Three or four blocks, or B=319, are clearly worse.
- Every from-scratch fit lands in a much worse basin than experiment 3's.

**Calibration: measured.** Relaxing experiment 3's 8-bit R=1 solution to float
(`exp3float.npz`) gives an error of 0.0113. Experiment 3's basin is therefore about twice as
accurate as anything I fitted from scratch, so effort and initialization were the bottleneck.

**Splitting experiment 3's solution: measured, failed.** Co-clustering its vectors and
columns into two blocks (`split_blocks.py`) gave an error near 1.5 even at full rank.

**Pruning from experiment 3's basin: measured.** `prune_blocks.py` starts both blocks with
the full 719-vector dictionary and co-clustered column halves. It then drops vectors, with a
variable-projection Adam refit between drops. Float error₂ by block size:

| Log | Schedule | Criterion | B=647 | 551 | 503 | 479 | 455 | 447 | 415 | 383 |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `pr_r2` | drop 24, 300 steps | diag + row norm of K | 0.0107 | 0.0153 | 0.0196 | 0.0225 | 0.0268 | stopped | | |
| `pr2_r2` | drop 8, 600 steps | diag + row norm of K | 0.0096 | 0.0119 | 0.0139 | 0.0154 | 0.0170 | 0.0176 | stopped | |
| `pr3_r2` | drop 8, 600 steps | leverage | 0.0090 | 0.0104 | 0.0119 | 0.0130 | 0.0146 | 0.0148 | 0.0171 | 0.0201 |

The leverage criterion is 2(NSN)_bb/Q_bb − (N_bb/Q_bb)². The extra capacity of two full
dictionaries pushes the error below the starting 0.0113 at first (0.0089 at B=695). Below
about B=550 the error rises steadily.

**8-bit quantization: measured, failed the 20 mHa cap.** The quantization pipeline is
`quanttest.sh`: `blockanneal.py` (continuous to lattice angle annealing), then
`block_coordinate.py` (integer descent with exact block weights, curvature-aware scoring and
pair moves).

| Input | Float error | After anneal | After descent | Chemistry gain if it passed (estimate) |
| --- | ---: | ---: | ---: | ---: |
| B=479 snapshot | 0.0129 | 0.118 | **0.0250** | 1.14× |
| B=447 snapshot | 0.0148 | 0.128 | ~0.0278 (stopped at step 1,755) | 1.20× |
| Experiment 3's basin, R=1 (calibration) | 0.0113 | 0.090 | 0.0227 | n/a (experiment 3's own 8-bit result is 0.0196) |

Integer descent took several rounds of tool development:
- Gradient-sign moves stalled at 0.163 (`qc_test`).
- Curvature-aware scoring reached about 0.055 (`qc2_test`).
- Adding pair moves reached about 0.047 on an R=2 fit (`qA`, `qA2`).
- An R=1 control stalled the same way (`qR1`, 0.063).

**Conclusion.** My quantizer is about 0.003 weaker than experiment 3's on the same basin.
Even allowing for that, R=2 at B ≤ 479 lands around 0.022 mHa, over the cap. Passing seems to
need B > ~500, worth about 1.1× at most. Experiment 3 independently recorded block fits
around 31 mHa after 8-bit rounding. This rules out *my* fits and quantizer, not the
architecture: a better quantizer or a much better float basin could change it.

### Other chemistry levers

- **Angle precision (estimate).** The cost model picks its QROAM block size k as the largest
  k with temporary qubits below brot·(N − 1).
  - At brot 9 that gives k = 5 and Q ≈ 823–856, worse for every shape scanned.
  - brot 7 would cut the angle circuits (K = 5 at J = 15) and the Givens rotations, about
    1.07× on T. I did not build or quantize a 7-bit fit; the expected doubling of rounding
    error is an argument, not a measurement.
- **Smaller single-block shapes (estimate plus experiment 3's evidence).** R=1 with B ≤ 511
  would be worth 1.38–1.49× in the cost model. I did not attempt such a fit. Experiment 3
  needed 719 vectors to pass at 8 bits, its solution sits at 19.995 mHa, and removing one
  column took it to 20.0095 mHa.
- **λ_eff (measured decomposition).** For the current solution,
  2Λ − E_gap = 559.13 (fixed by the verifier's shift) + 2Σλ⁺ (≈3) + a mixing term (48.45).
  The floor is λ_eff ≥ 46.79, against 48.895 now, so at most 1.045×. Experiment 3 already ran
  orthogonal column rotations, which change only the mixing term, and mixing-penalty refits.
  `constrained_block.py`, the block version of experiment 3's λ refit, was validated only by
  reproducing about 48.889 on experiment 3's solution.
- **Angle-loading circuits (analysis, wrong).** I concluded the 128-qubit per-oracle cap
  pins them at about 2·K·rows Toffolis, with J ≤ 14 and K = 6. Codex later beat this; see
  [Corrections](#corrections-in-light-of-codexs-results).

## DQI

All DQI work here was profiling and analysis of experiment 3's syndrome-consuming
Berlekamp–Massey (BM) decoder. No new decoder was built.

- **Reproduction (measured).** `dqi/cand/base.kmb` reproduces experiment 3's circuit:
  1,429 qubits and 645,755 Toffolis.
- **Peak profile (measured, instrumented `dqi/stream_capture.py`).** The peak of 1,428
  qubits is at BM round 68. Per-round peaks rise steadily, from 1,002 at round 0 to 1,428,
  as the locator fills in. At the peak the live qubits are:

  | Item | Qubits |
  | --- | ---: |
  | Syndrome register | 840 |
  | Locator σ | 432 |
  | Bank | 36 |
  | Stack | 19 |
  | Held registers | 23 |
  | Exception handling (active, filled, stop index 9, park index 6, park 12) | 29 |
  | Round temporaries | about 49 |

- **BM Toffolis (measured).** The forward pass costs 159,664 and the reverse 133,839, for
  293,503 in total.
  - Even rounds reverse at 0.65× their forward cost; odd rounds at about 1.02×.
  - Odd rounds absorb syndrome S_l in place (`mon_slot` plus XORed products). The reverse
    must regenerate S_l = σ_l⁻¹(Δ + R), which needs the same products. The apparent
    31K-Toffoli gap is the price of the qubit savings, not a fixable inefficiency
    (analysis).
- **Forney sweep (measured).** The per-point kernel costs 80 Toffolis: zero test 11, tower
  inverse-parity 66, final parity 3 (`recursive_forney.py` self-test).
- **Error-weight distribution (measured from the verifier's formula).**
  - Weight 35 occurs in 2.77% of shots and weight ≥ 34 in 12%, so capping the decoder below
    35 errors is ruled out.
  - The repetition factor is 1.191 at zero failures and 1.223 at the current 8.
- **Estimates and analysis (no circuits built):**
  - A Berlekamp–Massey state needs about 72 + L symbol slots, roughly 1,300 qubits.
  - Removing B or the discrepancy window needs a degree-35 Euclid, about 2,450
    multiplications.
  - An in-place elimination decoder needs about 57K multiplications per direction.
  - A qd-style in-place variant fails about 26% of the time over GF(4096).
  - Upstream's Ψ = 1/σ′ sweep (about 23 Toffolis per point) needs 35 more symbol slots
    while B is live, so the qubits² penalty outweighs the T gain.
  - Interpolation could move 70 sweep points onto free Z gates, about 0.9% of DQI T. Not
    implemented.
- **Literature.**
  - Gu et al. ([arXiv:2610.22367](https://arxiv.org/abs/2610.22367)): their best
    Reed–Solomon decoder has Q²T = 1.48e12, against 1.32e12 for experiment 3's decoder.
  - The competition paper ([arXiv:2610.23832](https://arxiv.org/abs/2610.23832)) withholds
    its best results.

## Inverse square root

No circuits were built beyond reproducing experiment 3's (`inv/base.kmb`, with stage marks
added in `inv/binary_centered_invsqrt.py`).
- **Measured.** The stages cost: normalization 210, lookup 308, Horner 1,043, output
  shift 228 Toffolis. The circuit also has 34 rotations (about 250 Toffoli-equivalents) and
  uses 98% of the error tolerance. Setting `qubit_limit` below 94 fails to build in that
  design.
- **Estimates:**
  - A cubic Horner would cost about 1,760 Toffolis, against 1,043 plus a 308-Toffoli lookup
    today.
  - A linear fit needs about 4,096 table entries.
  - Phasing each term separately needs hundreds of extra rotations.
  - Exponent-dependent output rotations multiply the rotation count.
  - Measurement-based cleanup of the output still needs its bits to correct the phase.

## Corrections in light of Codex's results

I was wrong about the following, and told the user these things as conclusions. The round
report (`../REPORT.md`) records Codex's natively verified results, which disproved them.
[MISTAKES.md](MISTAKES.md) has my original wording and the faulty assumption behind each.

1. **I said the chemistry angle loading was at its floor. It wasn't.** Codex's signed angle
   gauges and eleven 7-word oracles cut Q from 526 to 478, a 6.42% gain on chemistry.
2. **I said invsqrt was near optimal, with 94 qubits a hard limit. It wasn't.** Codex reached
   90 qubits and a 4.07% gain with affine lookup compression and delayed coefficient loading.
3. **I said DQI was near a qubit floor of about 1,300, and that compact decoders could not
   come out ahead. Both were wrong.** Codex's compact Euclidean decoder peaks at 961 qubits
   and nets +0.90%. My Toffoli estimate for compact decoding was several times too high.
4. **I said the ceiling was about 9.05, and that 9.4 was unreachable incrementally.** The
   first claim was wrong: Codex reached 9.3677 incrementally. Before that I had told the user
   "9.5–10 is realistic", which led them to set the 9.4 goal.
5. **I said chemistry λ_eff was "exhausted".** Codex lowered it from 48.895 to 48.594.

Codex's results agreed with mine on the λ_eff bound (46.786 against my 46.79), on smaller or
blocked dictionaries failing quantization, and on the cubic invsqrt being worse.

## What would plausibly move the score

- **Chemistry.** Port Codex's signed-gauge batching, then return to smaller or blocked
  shapes. The quantizer is the weak link. The best basin anyone has is experiment 3's (0.0113
  float); work from it rather than from scratch.
- **DQI.** Look for decoders that trade width for T without paying a large repetition
  factor. Codex's unfinished Wronskian padding adapter (28,069 → 7,140 Toffolis per
  direction, components only) is the open lead.
- **Invsqrt.** Start from Codex's 90-qubit circuit, not from experiment 3's 94-qubit one.
