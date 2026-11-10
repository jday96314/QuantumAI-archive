# Round 5 chemistry quantization-aware factorization

**Final selection: no replacement. The preferred 567M and stretch 506M targets
were not reached.** A smaller R1/B700/C360 eight-bit factorization passes full
native verification, but costs **825,280,181.6200**, so it is rejected. The
best new fully verified candidate costs **728,405,368.7660**, also worse than
the preserved round-4 **708,547,388.9553** incumbent. Neither is a task-cost win.

Work is confined to `initial_experiments/5/chem-qfit`. Earlier rounds and host
assets are read-only. No other round-5 directory was consulted. The selected
combined reference remains 9.3677221531 (inverse 182956.5272, DQI
1598204229385, chemistry 708547388.9553). No replacement is selected.

## Calibration

The actual round-4 Codex-selected chemistry is
`initial_experiments/4/codex/candidates/norm_cycle2_signed`, as bound by
`submission.summary.json`, not the older signed-gauge families mentioned in
stale notes. `results/calibration.json` preserves all source SHA256s.

Native physical metrics reproduce exactly: error1 0.011986315442796802,
error2 0.0199953044905532, Lambda 303.5283147431962, gap
3.915157472365536, lambda_eff 48.59417497035283. Signed decode/re-encode
changes **zero words** and zero projector entries. Metrics are identical.
Exact pair-space error is 0.019995303803338634 (native difference 6.9e-10).
The reference spectrum contains negative eigenvalues (min -0.00029765);
new optimizers use the full signed spectrum rather than clip it to PSD.

Fixed incumbent words, unregularized rank-350 weight projection: error2
19.5330304 mHa; ridge 1e-5: 19.5345981 mHa. Thus the ridge alone does not
explain the earlier failed control. Claude's `L_r1b719c350` float basin has
21.0283214 mHa true error before quantization, and its `qd_r1cal` integer
result is 22.6571803 mHa. Neither is the accepted incumbent basin.

Matched starting basin diagnostics (mHa; full signed target; four exact block
weight sweeps; these are numerical metrics, not native circuit passes):

| Shape | Float/saved | Direct rounding | Exact integer-word weight refit |
|---|---:|---:|---:|
| R1/B719/C350 Claude L | 21.0283 | 5980.17 | 142.048 |
| R2/B479/C175 snap | 13.1453 | 210.991 | 145.433 |
| R2/B447/C175 snap | 14.8921 | 252.602 | 170.660 |
| R2/B351/C175 L | 21.1151 | 553.415 | 180.195 |
| R1/B719/C350 Claude integer control | 22.6572 | 22.6572 | 22.6558 |

`B` is per block: R2/B479 means 958 two-body vectors. R2/B351 is excluded
from discrete experiments unless a float basin with adequate margin appears.

Reproduce from this directory:

```bash
source env.sh
"$QFIT_PY" calibrate.py > results/calibration.log 2>&1
"$QFIT_PY" "$QFIT_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  "$QFIT_ROOT/initial_experiments/4/codex/candidates/norm_cycle2_signed" \
  --salt round5-chem-qfit-calibration > results/incumbent_native.json
```

`core.py` uses the native Givens convention, the prior signed encoder, full
pair-space tensor residual, exact QR block PSD/rank-C projections, and retains
warm weights if an update would increase the true residual. Imported host
DFTHC code is Apache-2.0; prior Codex/Claude/AGY work is attributed by path.

## Discrete implementations and controls

1. `freeze.py`: tensor-error feedback selects vectors or individual angles to
   round and freeze, then retrains the remaining continuous coordinates and
   weights on the full tensor residual. Exact float64 word decoding and block
   PSD/rank-C residual refits evaluate every checkpoint. Tested least-impact
   first, random angle order, and hardest-vector first. Initial whole-vector
   freezing was too coarse. On the identical `snap_B447.npz` float basin, the
   initial least-impact run gave **116.3245 mHa** and the corrected hardest-first
   run **27.5371 mHa**, with different retraining budgets and gauge-continuity
   handling. This is not an isolated ordering comparison; the improved run still
   fails. Random angle freezing on `snap_B479` reached **40.5166 mHa**, versus
   **101.1805** for coarse vector freezing and **107.8607** for least-impact
   angle freezing. None beats the previous best B479 quantizer near 25 mHa.
2. `discrete.py`: exact integer moves couple 1/2/4/8 words within a vector and
   combine moves across vectors. Batched nonlinear decoded-feature differences
   screen moves against the true tensor residual and dictionary tangent space;
   exact rank-C weight refits decide acceptance. It reproduces the accepted
   basin and reaches **19.5163741 mHa** in 40 iterations. This is a numerical
   improvement of the control, not a new scored submission. On B700 it was
   slower and less effective than the existing curvature solver. On Claude's
   identical B479 integer basin, a fixed-word 100-sweep refit control gives
   **24.9572812 mHa**; new moved words under the same refit give
   **24.9499002 mHa**. That small difference does not solve quantization.
3. `exchange.py`: constructs integer projector vectors from leading tensor
   residual modes, screens whole-row exchanges by their exact common-span
   projection-energy change, then accepts only after true rank-C tensor refits.
   It improves pruned B671 from **24.5053 to 23.7993 mHa**; none of the tested
   exchanges helped B700 or Claude's B479 basin. This is a finite proposal
   screen, not a proof that no beneficial exchange exists.

`legacy_block_exact.py` and `anneal_exact.py` adapt the prior block coordinate
and coupled QUBO controls, with the full signed target. `legacy_block_control`
runs the unmodified round-4 Claude optimizer from the actual accepted words:
its logged PSD-target result is 19.5239046 mHa; recomputation against the full
signed target gives **19.5262560 mHa**. Thus its known-basin control does pass.
The 300-step round-3 STE control starts at 19.9953 mHa, jumps above 30 mHa after
enabling angle updates, and keeps its initial snapshot as best. No scalar angle
rounding loss was used as an acceptance test.

`prepare.py`, `greedy_prune.py`, and `split_integer.py` generate smaller exact
word dictionaries from the accepted basin. Sequential pruning recomputes
leverage after deletion. Integer column splitting avoids the large initial
rounding shock of the two-block float fits; its B479 integer descent reached
about **22.62 mHa**, still invalid. See `results/numerical_summary.json` for
recomputed errors of every main saved numerical candidate, and
`results/stopped_runs.json` for intentionally truncated plateaued runs. Their
saved K checkpoints were converted to explicit weights before handoff.

Two corrections are preserved in the implementation: exact eigensolves align
their column gauge to the warm weights by orthogonal Procrustes, avoiding
unnecessary optimizer-state and normalization changes; prior unfinished
checkpoints may store K rather than w, with either one or multiple blocks.
The loader reconstructs the positive numerical rank explicitly. The native
float32 shifts are always obtained through `load_inputs`; they are never fitted.

A final matched ablation isolates gauge continuity on `snap_B447`, keeping the
same seed, warmup=300, 300 retraining steps per batch, batch=40, hardest-first
order and all other settings: aligned projections finish at **29.5274297 mHa**;
unaligned projections finish at **29.1939450 mHa**. Alignment preserves a useful
column convention but does **not** improve quantization in this control. It
cannot explain the large earlier accuracy loss. Both exact numerical artifacts
and logs are retained as `ablation_align_B447` / `ablation_noalign_B447`.

## Host resource screens

`screen_costs.py` builds every full signed-angle oracle in each screen; the
following scores use actual host state-prep arithmetic and complete circuit
static Toffolis, but **assume the incumbent lambda_eff**. They are resource
screens, not native verified candidate scores. All use brot=8, J=7, chi=10,
and the host chooses k2=4. The full JSON includes inner temporary peak, all
oracle widths/counts and the normalization needed to reach 567M.

| R | B per block | total R*B | C per block | T total | Q total | Cost at incumbent lambda (M) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 511 | 511 | 350 | 23433 | 460 | 523.805 |
| 1 | 639 | 639 | 350 | 27794 | 478 | 645.599 |
| 1 | 671 | 671 | 350 | 28878 | 478 | 670.779 |
| 1 | 700 | 700 | 350 | 29915 | 478 | 694.866 |
| 1 | 719 | 719 | 350 | 30504 | 478 | 708.547 |
| 2 | 447 | 894 | 175 | 26171 | 461 | 586.281 |
| 2 | 479 | 958 | 175 | 27981 | 461 | 626.828 |
| 2 | 351 | 702 | 175 | 21863 | 460 | 488.711 |

C160 and C190 were also screened for B447/B479; low C saves inner cost but
raises rank truncation error. R2/B351 had 21.1151 mHa in the available float
basin and was excluded from further discrete optimization. The B511 probe
did not establish adequate continuous margin, so its failed run is a numerical
initialization/training result, not evidence of a B511 capacity floor.

`oracles.py` preserves the selected XOR-difference interface and signed gauge.
At A=1034 (R2/B479), the 10 high address bits need 52 clean qubits for monomial
phase repair, exceeding the 49 difference payload bits. Padding with three
zero qubits yields **127Q**, within the cap. The unpadded builder's assertion
was a construction limitation, not a physical failure. Unsigned legacy words
must first be re-gauged by decoding and signed re-encoding; projectors are
preserved. `results/oracle_forward_reverse.json` records all-address tests with
random incoming words, correct values, restored addresses, clean scratch and
phase, in both directions: 1590 shots/direction/oracle for the incumbent and
2068 for the padded B479 family. This is an oracle-component check, not a
physical pass for B479.

## Physical checks and exact task costs

The generalized AGY offset solver in `physical.py` uses native normalized X,
native fixed shifts, all R*C offsets, and the actual SOS gap. Its gap/Lambda
match the incumbent native metrics within 4e-10; four directional gradient
checks have maximum error 7.4e-9 (`results/offset_semantics.json`). It is a
convex gap optimizer, not a global task-cost solver. A below-floor minimizer
must be moved to a strict accepted gap. One-body words are optimized against
the actual weighted one-body matrix, not angles, and then signed re-gauged.

Orthogonal column-gauge rotations use the prior read-only Codex `jacobi_l1`
binary and preserve the combined weight/identity Gram matrix; recorded errors
are below 1e-9. `constrain.py` adapts the prior constrained weight/offset refit
only for a smaller candidate, using the physical error and SOS gates. It did
not produce a winning cost. No normalization-only incumbent is proposed.

| Candidate | Native status / full shots | error1 (mHa) | error2 (mHa) | lambda_eff | T | Q | Native cost |
|---|---|---:|---:|---:|---:|---:|---:|
| Round-4 selected | PASS / 50880 | 11.986315 | 19.995304 | 48.594175 | 30504 | 478 | 708547388.9553 |
| Greedy B711/C350 | PASS / 50368 | 14.937006 | 19.828293 | 56.117926 | 30408 | 478 | 815675394.3637 |
| Integer B700/C360 | PASS / 49664 | 14.941592 | 19.842585 | 56.853516 | 30368 | 478 | 825280181.6200 |
| Constrained B711/C350 | PASS / 50368 | 14.999636 | 19.990336 | 50.113806 | 30408 | 478 | 728405368.7660 |

Each PASS is the unchanged native verifier with a distinct explicit validation
salt, both physical residual bounds, strict Lambda and gap floors, the full
forward/reverse ANGLE family and oracle-width cap. JSONs are respectively
`incumbent_native.json`, `greedy711_native.json`, `b700c360_native.json`, and
`constrained711_native.json`. Every new score is worse and is rejected.
Complete NPZ/circuit directories are `candidates/greedy_B711/`,
`candidates/legacy_exact_B700_C360/`, and `candidates/constrained_B711/`.
Additional B700/C365 and C375 builds are resource/physical screens only; they
were not native verified and are not selected.

For B700/C360 the inner preparation is **17414 T and 337 temporary qubits**,
versus the incumbent's 17333 T/337 temporary qubits. The smaller angle table
saves gates, but extra C and worse physical normalization erase the gain.

## Failure classification

* **Encoding/calibration:** exact round trip and native incumbent reproduction
  pass. The accepted and Claude default trees are identical. The ridge changes
  the incumbent refit by only 0.0016 mHa. The failed R1 calibration used a
  different already-invalid float basin; the matched accepted control passes.
* **Identical float basins:** `snap_B479` and `snap_B447` have 13.1453 and
  14.8921 mHa full-target errors; rounding/refitting gives 145.433/170.660 mHa.
  The large deterioration is quantization and subsequent discrete optimization,
  not a lack of continuous shape capacity. Matched single-block float/rounded/
  refitted probes are in `results/matched_B*.log` and saved warmup NPZs.
* **Fixed dictionary obstruction:** the particular legacy quantized B447 words
  have a **20.5962034 mHa** relaxed joint-span lower bound. Even unlimited
  shared real coefficients cannot rescue that dictionary by weights alone.
  This does not exclude different integer words for the same shape.
* **Unresolved block/rank optimizer:** B479's joint-span bound is
  **17.8537093 mHa**, while 12 block PSD refit sweeps give 24.96044 mHa at
  C175, 24.01126 at C190, and 22.93562 at C225. Increasing C alone does not
  provide a passing candidate. The joint-span relaxation permits cross-block
  coefficients absent from the submitted model, so its bound is not an
  achievable block-model claim. These experiments do not establish a general
  quantization floor or globally optimal block fit.
* **Cost failure:** smaller fully passing models exist (above), but fail the
  actual lambda_eff*T_total*Q_total objective. They are not promoted, and the
  combined reference stays 9.3677221531.

Matched single-block probes use 3000 Adam warmup steps, exact float64 weight
projection, and the same trained weights for direct rounding (mHa):

| Shape | Continuous | Direct rounded | Exact integer-word refit |
|---|---:|---:|---:|
| R1/B511/C350 | 48.0525 | 229.4015 | 165.3178 |
| R1/B639/C350 | 15.1118 | 145.8143 | 101.5521 |
| R1/B671/C350 | 13.5986 | 122.9565 | 90.1967 |
| R1/B700/C350 | 13.3231 | 119.3779 | 87.3839 |

These show why continued float-only fitting is not the answer for B639–700.
Direct integer continuation from the accepted basin gives much better results
than re-quantizing these float basins, although it still misses the useful
smaller-shape thresholds in the tested searches.

## Reproduction of the main experiments

Run in this directory after `source env.sh`. Scripts write only here; imported
earlier-round modules/binaries are read-only. Numerical searches use the local
CUDA torch environment; no host scoring code is changed.

```bash
"$QFIT_PY" prepare.py > results/prepare.log 2>&1
"$QFIT_PY" screen_costs.py > results/shape_costs.log 2>&1
"$QFIT_PY" discrete.py --init candidates/incumbent_control.npz \
  --out candidates/discrete_control.npz --steps 40 > results/discrete_control.log
"$QFIT_PY" freeze.py \
  --init "$QFIT_ROOT/initial_experiments/4/claude/chem/snap_B447.npz" \
  --out candidates/hardfreeze_B447.npz --warmup 300 --steps 500 \
  --batch 24 --order hardest --lr .03 > results/hardfreeze_B447.log
"$QFIT_PY" rank_diagnostics.py > results/rank_diagnostics.log
"$QFIT_PY" check_oracles.py > results/oracle_forward_reverse.log
"$QFIT_PY" physical.py candidates/legacy_exact_B700_C360.npz \
  --align --jacobi 12 --build > results/physical_rank_aligned.log
"$QFIT_PY" "$QFIT_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  candidates/legacy_exact_B700_C360 --salt round5-chem-qfit-b700c360 \
  > results/b700c360_native.json
"$QFIT_PY" physical.py candidates/constrained_B711.npz \
  --preserve-offset --jacobi 12 --build > results/physical_constrained.log
"$QFIT_PY" "$QFIT_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  candidates/constrained_B711 --salt round5-chem-qfit-constrained711 \
  > results/constrained711_native.json
"$QFIT_PY" summarize.py > results/numerical_summary.log
```

The original submitted candidate NPZs and new exact-word search checkpoints
are retained, so native verification does not require rerunning numerical
search. A fresh stochastic/numerical search may differ in final low-order
floating-point digits; native JSONs bind the exact delivered artifacts by
submission digest. See `artifact_manifest.json` and `selection.json`.

Final audit: all four native digests match the delivered circuit/NPZ bytes.
Rebuilding the incumbent through `oracles.py` produces all eleven circuits and
the NPZ byte-for-byte identical to the preserved source; its archive SHA256 is
unchanged. Evidence: `results/audit.log`, `results/incumbent_rebuild.json`, and
`artifact_manifest.json`. All experiment processes are terminal; intentionally
stopped numerical runs retain their exact-word checkpoints and stop reasons.
