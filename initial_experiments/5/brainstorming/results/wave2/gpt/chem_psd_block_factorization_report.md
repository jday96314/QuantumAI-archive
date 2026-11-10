# Round 5 continued: block-PSD chemistry factorization

**Complete; no replacement. The 1.3x target was not reached.** The selected round-4 chemistry remains R1/B719/C350, eight bits, T=30,504, Q=478, lambda_eff=48.59417497035283, cost **708,547,388.9553173**. Fresh unchanged-host verification passes 50,880 shots with an independent salt. No numerical search or experiment process remains running.

The main positive result is a continuous R2/B351 factorization at **C189, 19.978781 mHa**, with conditional T/Q **22,516/461** (1.4047x at incumbent normalization). This is not an integer or scored submission. The exact rounded dictionary has a **140.196114 mHa joint-span lower bound** and a roughly 145.48 mHa block-PSD optimum. Increased columns help the continuous representation but cannot rescue those decoded vectors.

The best newly constructed exact eight-bit dictionary is the residual-augmented **R2/B375** pivot: full-PSD error **51.459643 mHa**, explicit dual lower bound above 51 mHa. Its actual rank-C192 native word screen is **53.181669 mHa**, also invalid. These are fixed-dictionary failures, not global limits on these shapes or on chemistry factorization.

## Scope and calibration

All writes are confined to this directory. Inputs are the unchanged host, the round-4 accepted candidate, completed `chem-blocks/RESULTS.md`, `chem-blocks/jad_search/RESULTS.md`, both resource-frontier implementations, and completed `chem-qfit/RESULTS.md` and its referenced round-4 checkpoints. No concurrently active round-5 agent directory was used.

`core.py` constructs the reference directly through native `load_inputs`: float32 FCIDUMP integrals, native float32 fixed S, and fixed beta_1. It uses `P[b,pq]=u[b,p]*u[b,q]*fac[pq]`, with fac=sqrt(2) off diagonal. The signed target is never clipped to PSD. Its negative spectral component is 0.297649784 mHa. Direct four-index, pair-space and native incumbent errors agree to about 2e-16 Ha.

`results/calibration.json` records the native control and the analytic single-block PSD control; `envelope_gradient_control.json`, `hessian_control.json`, `offset_controls.json` record finite-difference checks. The final Hessian implementation uses symmetric-pair reshuffling without a dense N^4 intermediate. Its directional error is about 8.7e-12. Native gap/Lambda formulas match the incumbent within 1e-8; fixed shifts are never fitted.

## Fixed-dictionary convex optimization

For each block, QR of P_r.T gives an orthonormal row dictionary Q_r and triangular R_r. Optimize S_r>=0 in `Hhat=sum Q_r.T S_r Q_r`; atom covariance is `K_r=R_r^-1 S_r R_r^-T`. This avoids projecting PSD in an incorrect atom metric. The solver caches `Q_r Q_s.T` and `Q_r H Q_r.T`. FISTA uses the safe Lipschitz bound R, adaptive restart, and full-capacity eigenvalue projection. Rank-C comparisons use exact cyclic block minimization by positive top-C spectral projection.

Loss contractions are checked against direct pair-space tensors and explicit factor weights. All generator modes can activate; extra zero factor columns do not remain artificially inactive. Rank-constrained fits are nonconvex feasible upper bounds, with multiple warm starts and rank-aware projected-gradient diagnostics. They are not global rank-optimality claims.

For E=Hhat-H, the dual tensor `D=E+delta*Pi_union` is feasible when every `Q_r D Q_r.T>=0`. It gives squared-error lower bound `-||D||^2-2<D,H>`. During optimization a conservative cached formula uses an upper bound on union rank. `certify.py` independently reconstructs dense final dual tensors, checks positive feasibility slack and direct weight reconstruction. Certificates use floating-point safety slack, not interval arithmetic. The reported excluded margins are far larger than numerical reconstruction errors.

### Experiment A: prior dictionaries

| Dictionary | C88 | C128 | C160 | C175 | C192 | C224 | Full PSD feasible / dual lower (mHa) |
|---|---:|---:|---:|---:|---:|---:|---:|
| prior_r4b180 | 54.677 | 50.736 | 50.570 | 50.570 | — | — | 50.570 / 50.462 |
| prior_r4b224 | 31.725 | 28.032 | 27.579 | 27.564 | — | — | 27.564 / 27.328 |
| prior_r2b351 | 161.346 | 55.338 | 32.005 | 27.263 | 25.671 | 24.422 | 23.898 / 23.799 |

Extra C cannot rescue any of these exact prior dictionaries. The common-span relaxation further distinguishes their obstructions:

| Dictionary | Unrestricted joint-span PSD lower bound (mHa) | Interpretation |
|---|---:|---|
| prior_r4b180_psd | 33.015 | Union span already excludes 20 mHa |
| prior_r4b224_psd | 16.125 | Joint span permits <20; within-block PSD geometry excludes it |
| prior_r2b351_psd | 18.807 | Joint span permits <20; within-block PSD geometry excludes it |
| claude_r2b351_psd | 14.559 | Continuous full block-PSD fit passes |
| claude_r2b351_q8_psd | 140.196 | Quantized union span itself excludes 20 mHa |
| integer_newton_r2b351 | 47.250 | Integer moves improve the span but do not pass |
| residual_r2b375_psd | 45.060 | Residual additions improve it further but do not pass |

The joint-span relaxation permits cross-block coefficients and therefore is only a lower bound, not an achievable block representation. SVD dimensions and smallest retained singular values are recorded in `span_diagnostics.json`. Native unrestricted PSD rank must be at least 290; the resource CSV explicitly excludes shapes failing this necessary test.

### Stronger continuous basin and column threshold

The completed qfit report identifies `4/claude/chem/L_r2b351c175.npz`, whose continuous error is about 21.115 mHa. Its fixed-dictionary full-PSD interval is approximately **18.428–18.538 mHa**. Cold truncation of the relaxed covariance sometimes underperformed the original rank-C weights; both starts were retained and refined.

| C | Best demonstrated continuous error (mHa) |
|---:|---:|
| 160 | 26.844676 |
| 175 | 21.114437 |
| 184 | 20.297812 |
| 188 | 20.037224 |
| 189 | 19.978781 |
| 190 | 19.922480 |
| 191 | 19.867131 |
| 192 | 19.816275 |
| 200 | 19.465370 |
| 208 | 19.208374 |
| 224 | 18.881694 |

**C189 is the smallest demonstrated passing continuous rank.** C188 reached 20.037224 mHa in the tested refinement; that is not a certificate that C188 or a smaller C cannot pass with another rank optimizer. No physically viable C exists for the exact decoded dictionary, because its full-rank relaxation and common-span bound both fail. `claude_r2b351_*_capacity.json` preserves complete curves and warm/cold evidence.

![Rank-versus-error curves](results/rank_curves.png)

## Experiment B: different geometry and exact lattice moves

Two new non-JAD initializations were tested: (1) an accepted-covariance graph with high-coupling vectors shared across ensembles; (2) individual target-generator eigenvectors, spatial-participation grouping, and nonorthogonal simultaneous projector pivots. The graph R4/B224 start shares 59 vectors and covers all 719 incumbent atoms in its union. It nevertheless has a very poor block-PSD optimum, showing that coverage alone is insufficient.

`geometry_trust.py` improves unit vectors by PSD variable projection with curvature history and a physical preconditioner. A feasible-covariance Armijo test guarantees physical descent before the next covariance solve. An ordinary inexact-inner L-BFGS control terminated prematurely and is retained as a numerical optimization failure. The trust searches use about 100 accepted geometry steps, thousands of convex covariance updates, and independent long final PSD solves. GPU warm runs were continued on CPU after a timing control; exact saved checkpoints bind the reported results.

| Geometry | Initial full-PSD error (mHa) | Final feasible / dual lower (mHa) | Accepted steps |
|---|---:|---:|---:|
| trust_r4b224_cpu | 26.240 | 23.865 / 23.843 | 102 |
| trust_prior_r2b351_cpu | 23.898 | 22.137 / 22.082 | 102 |
| trust_eigen_r2b351 | 542.220 | 38.518 / 38.478 | 102 |
| trust_graph_r4b224 | 1578.125 | 168.094 / 168.094 | 103 |

The first row continues the preserved GPU warm geometry (the original R4/B224 control was 27.564 mHa). Geometry budgets and exact stop reasons are in `geometry_budgets.json`. These terminal fixed-dictionary certificates do not certify global geometry minima.

`integer_trust.py` then changes coupled Givens words across whole vectors, screens exact row tensor increments, tests combined decoded moves, and refits unrestricted PSD covariances. A second run uses validated per-vector physical Newton curvature. A bounded projected-trial lookahead pilot is also retained; its sparse logged history is explicitly marked as such. Every accepted move uses the true signed physical residual.

| Exact eight-bit geometry search | Initial / final full-PSD error (mHa) |
|---|---:|
| integer_claude_r2b351 | 145.476 / 57.344 |
| integer_newton_r2b351 | 68.084 / 54.718 |
| integer_projected_r2b351 | 68.084 / 64.583 |

One additional resource-screened family, R2/B375, appends 24 exactly encoded atoms per block to an integer checkpoint. Four additions are shared high-leverage residual vectors and twenty are block-specific. The pool combines incumbent vectors and individual eigenvectors of positive residual modes; the actual optimization still retains the full signed target. Adding these decoded nonorthogonal vectors lowers full-PSD error from 55.831479 to 51.459643 mHa. A convincing dual bound still excludes 20 mHa for this exact augmented dictionary.

![Geometry convergence](results/geometry_curves.png)

## Experiment C: decoded physical and normalization checks

Each row below restores rank C, decodes exact signed eight-bit words, refits weights, optimizes the actual SOS gap with all identity offsets, preserves native symmetry shifts, and evaluates native physics. Selected cases also receive column-gauge L1 optimization that preserves the combined weight/identity Gram matrix. Since every two-body dictionary fails even its relaxed PSD gate, further one-body polishing and constrained task-cost optimization cannot rescue these exact geometries. One-body results below are measured screens, not one-body impossibility claims.

| Screen | C | Continuous / decoded error2 (mHa) | Native error1 (mHa) | Lambda | Gap (Ha) | lambda_eff |
|---|---:|---:|---:|---:|---:|---:|
| claude_r2b351_c200 | 200 | 19.543 / 151.266 | 56.031 | 1373.088 | 5.265272 | 120.132 |
| integer_best_c200 | 200 | 55.812 / 55.795 | 51.376 | 1245.548 | 4.110607 | 101.109 |
| prior_r4b224_c128 | 128 | 28.032 / 140.519 | 61.962 | 1111.907 | 5.474636 | 110.203 |
| residual_r2b375_c192 | 192 | 53.202 / 53.182 | 50.096 | 1130.054 | 4.106325 | 96.249 |
| trust_eigen_r2b351 | 192 | 42.935 / 160.941 | 61.225 | 1708.612 | 5.601635 | 138.241 |
| trust_graph_r4b224 | 128 | 168.181 / 427.982 | 78.189 | 5075.510 | 13.147575 | 365.087 |
| trust_prior_r2b351_cpu | 192 | 26.183 / 156.501 | 55.101 | 1603.770 | 5.438466 | 131.964 |
| trust_r4b224_cpu | 128 | 26.342 / 134.403 | 56.091 | 1126.163 | 5.412962 | 110.283 |

Every new row fails both error gates; all pass the strict raw-Lambda/gap floors and gap-domain test. Independent decoded full-PSD refits in `decoded_psd_summary.json` have lower bounds above 20 mHa for all eight screened dictionaries. Rank-C physical fits can be worse than those relaxations because their optimizer is nonconvex. This difference is retained rather than mislabelled as a quantization lower bound.

## Complete host resource accounting

`resources.csv` contains 436 shape/precision points, with unmodified host chi10, automatic preparation workspace choices, explicit signed-ANGLE bounds and necessary signed-spectrum rank screens. Host preparation/Givens arithmetic is exact. ANGLE bounds are constructive worst cases; synthetic dense tables attain them and pass arbitrary-target XOR/phase/scratch tests under two independent seeds.

| Shape (eight bits) | T | Q | Conditional improvement at incumbent lambda | lambda_eff ceiling for 1.3x |
|---|---:|---:|---:|---:|
| R4/B180/C88 | 18,391 | 443 | 1.7897x | 66.8985 |
| R4/B180/C128 | 20,457 | 445 | 1.6017x | 59.8720 |
| R4/B180/C160 | 22,084 | 445 | 1.4837x | 55.4610 |
| R4/B224/C128 | 24,017 | 445 | 1.3643x | 50.9973 |
| R2/B351/C175 | 21,863 | 460 | 1.4498x | 54.1948 |
| R2/B351/C189 | 22,516 | 461 | 1.4047x | 52.5090 |
| R2/B351/C200 | 23,028 | 461 | 1.3735x | 51.3415 |
| R2/B351/C224 | 24,160 | 463 | 1.3035x | 48.7245 |
| R2/B375/C192 | 23,975 | 461 | 1.3192x | 49.3135 |

`resource_actual.py` also compiles **all eleven actual ANGLE circuits for each of the eight physical screens**, retained under `resource_circuits/`, and applies the exact host persistent-plus-maximum-temporary resource formula. Actual compiled T/Q agree with the bounds for these tables. The best two new physical products are still above 1.06 billion at their measured normalization, and every model is physically invalid. These products are not accepted task scores. No full native PASS is claimed for an invalid family.

Only the incumbent has complete native chemistry verification: `results/incumbent_native.json`, canonical digest `ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7`, 50,880 shots. Its original twelve files remain read-only. There is no new viable candidate requiring a new accepted native submission.

## Artifacts, limits, and reproduction

`selection.json` records no replacement. `artifact_manifest.json` hashes all research code, exact NPZs, JSON evidence, plots, actual circuit families, source inputs and baseline files. `parameters/` preserves continuous covariances, explicit rank-C weights, integer words, decoded fits, offsets and normalization gauge inputs/outputs. CSV, PNG and PDF curves are standalone artifacts.

The numerical conclusions concern the saved dictionaries and finite geometry searches. C189 is an empirical constructive threshold, not a proven minimum rank. Floating trajectories can vary by device and scheduling; physical metrics and certificates are recomputed from the exact hashed checkpoints. The failed initializations, premature L-BFGS control, finite lookahead pilot and CPU restarts are preserved. None establishes a global multiblock, eight-bit, normalization or chemistry limit.

Run commands from this directory after `source env.sh`; set `PSD_DEVICE=cpu` for CPU reproduction.

```bash
$PSD_PY calibrate.py
$PSD_PY resources.py
$PSD_PY resource_controls.py
$PSD_PY experiment_a.py
$PSD_PY dictionaries.py --R 4 --B 224 --shared 59 --strategy graph --tag overlap_graph_r4b224
$PSD_PY dictionaries.py --R 2 --B 351 --shared 32 --strategy eigenatoms --tag eigen_r2b351
$PSD_PY solve.py --input parameters/eigen_r2b351.npz --out parameters/eigen_r2b351_psd.npz --steps 2400
$PSD_PY geometry_trust.py --input parameters/eigen_r2b351_psd.npz --tag trust_eigen_r2b351 --steps 100 --inner 30 --radius .05
$PSD_PY capacity_scan.py --input parameters/claude_r2b351_psd.npz --tag reproduce_rank_curve --ranks 175 184 188 189 192 200 224 --full-steps 0
$PSD_PY integer_trust.py --input parameters/claude_r2b351_q8_psd.npz --tag reproduce_integer --steps 100 --inner 30
$PSD_PY integer_trust.py --input parameters/integer_newton_seed.npz --tag reproduce_newton --steps 100 --inner 30 --newton
$PSD_PY residual_pivot.py
$PSD_PY physical.py --input parameters/residual_r2b375_psd.npz --tag reproduce_residual --rank 192 --steps 1000 --jacobi 8 --build
$PSD_PY quantized_certificates.py
$PSD_PY resource_actual.py
PYTHONPATH="$PSD_DIR/.plotdeps:$PYTHONPATH" $PSD_PY plot_results.py
```

Host DFTHC imports retain the Apache-2.0 license at `from_host/verifiers/LICENSE.txt`. Signed encoding/oracles reuse read-only round-4 `angle_gauge_probe`, `build_selectswap_chem` and `monomial_unlookup`, plus completed chem-blocks builders. Exact integer-row increments and offset/gauge semantics are attributed in the local source. The plotting dependencies were installed only into `.plotdeps/`, with licenses included there; no host or global environment was changed.
