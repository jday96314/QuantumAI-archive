# Experiment 4 (claude)

## Status

No improved submission. `submission.zip` and `submission.summary.json` are experiment 3's
verified components, with an overall score of **9.0096**: invsqrt 190,712, DQI 1.6127e12 and
chemistry 757,163,489. Their verify JSONs are copied into `results/`.

The round had two goals:
- 14.0, later cleared;
- 9.4, not reached.

The Codex workspace from the same round reached 9.3677 (see `../REPORT.md`).
- [FINDINGS.md](FINDINGS.md) has every experiment and its result, labelled measured,
  host-verified or estimate.
- [MISTAKES.md](MISTAKES.md) records the limits and estimates I stated that Codex's results
  disproved: the angle-loading floor, the invsqrt and DQI qubit limits, the ceiling of about
  9.05, the "9.5–10 realistic" estimate and "λ exhausted". Read it before trusting any
  "floor" claim in these docs.

## Environment

- `source env.sh` sets `PYTHONPATH` (kickmix, the host verifiers, this directory) and `$PY`.
- GPU fits use `/mnt/data01/miniconda3/envs/torch2.13/bin/python`.
- Files from experiments 1–3 are only read, through absolute paths; nothing there was
  modified.

## Top level

| File | Purpose |
| --- | --- |
| `costscan.py` | Chemistry cost model: `cost(R, B, C, brot, J)` returns (T, Q, heuristic breakdown) from the verifier's state-preparation model plus the angle-circuit model |
| `package.py` | Packages and scores `submission.zip` (copied from experiment 3) |

## `chem/` — chemistry

### Float fitting

| Tool | Purpose |
| --- | --- |
| `shapefit.py` | Adam fit of an (R, B, C) block DF-THC shape, with periodic exact block-weight projection |
| `lbfit.py` | L-BFGS version of the same |
| `blockvp.py` | Variable-projection Adam: block PSD rank-C weights re-solved each step; optional rank annealing (`--anneal`) |
| `vpfit.py` | Variable-projection L-BFGS (stalls; not useful) |
| `init_eig.py`, `init_isdf.py` | Eigenvector and greedy ISDF block initializations |
| `split_blocks.py` | Co-clustering split of an R=1 solution into blocks (fails) |
| `prune_blocks.py` | Continuation from R=1 to R blocks by gradual pruning; `--criterion diag` or `leverage` |
| `evalfit.py` | Evaluation helper for float fits |

### Quantization and refits

| Tool | Purpose |
| --- | --- |
| `blockanneal.py` | Continuous to brot-bit lattice angle annealing for blocks |
| `block_coordinate.py` | Integer angle descent: exact block weights, curvature-aware scoring, pair moves |
| `requant.py` | Re-quantization helper |
| `quanttest.sh`, `quanttest_c.sh` | Run anneal then descent on one input (the second takes the column count as an argument) |
| `constrained_block.py` | Augmented-Lagrangian λ_eff refit under the error and gap gates (lr ~2e-7) |

### Build and verify

| Tool | Purpose |
| --- | --- |
| `pipeline.sh` | Descent, offsets, λ refit, one-body words, circuits, host verify |
| `optimize_offsets.py`, `quantize_onebody.py`, `build_global_chem.py`, `compress_chem.py`, `chemistry_qrom.py` | Experiment 3/2 build steps, re-pointed at experiment 3's data |
| `queue1.sh`, `queue2.sh`, `snapper.sh` | Batch drivers for shape scans and prune snapshots |

### Key data files

- `exp3float.npz`: experiment 3's 8-bit R=1 solution relaxed to float (0.0113).
- `pr3_r2.npz`, `snap_B{479,447,415,383}.npz`: leverage-pruned R=2 states.
- `an_*` / `qd_*`: anneal and descent outputs (`b479`, `b447`, and `r1cal` for the R=1
  calibration).
- Logs carry the same stem as their outputs; FINDINGS.md tabulates their final values.

## `dqi/` — DQI (profiling only)

These are copies of experiment 3/2 DQI builders, with `import_published.py` re-pointed at
experiment 3's vendored upstream.
- `cand/base.kmb` reproduces experiment 3's decoder.
- `stream_capture.py` is instrumented. Its outputs in `parameters/` are:
  - `<tag>_phasepeaks.json`: per-phase peaks;
  - `<tag>_tofs.json`: per-phase Toffolis;
  - a printed peak snapshot of claim owners, plus the allocation phases of qubits with no
    owner.

  Run it with `$PY stream_capture.py snapA_offset3_select_noscale_compact`. The output tag
  selects the variant options.
- `model/bmid.py` checks Berlekamp–Massey identities classically.
- `recursive_forney.py` self-tests the 80-Toffoli per-point kernel.

## `inv/` — inverse square root (profiling only)

- `binary_centered_invsqrt.py` is experiment 3's builder with stage marks added.
- `build_inv.py` rebuilds `base.kmb`, the 94-qubit circuit.
