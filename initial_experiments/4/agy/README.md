# Experiment 4 (AGY): Quantum Circuit Compilation & Benchmark Optimization

This directory contains all experimental pipelines, automated compilation scripts, verification logs, mathematical proofs, and submission deliverables for the Verifiable Quantum Circuit Compilation Competition.

---

## 1. Primary Documentation & Reports

The markdown documentation in this directory comprehensively captures all experiments, mathematical analyses, and findings across all three benchmarks:

1. [**`EXPERIMENT_REPORT.md`**](EXPERIMENT_REPORT.md): **The Master Experiment Report.** Contains exhaustive documentation of all experiments across Task 1 (InvSqrt width sweeps $Q=86\dots 93$), Task 2 (DQI Fitzpatrick EEA sweeps, 23-Toffoli $\Psi$ sweep analysis, break-even width derivations, and in-place rational reconstruction architecture), Task 3 (Chemistry analytical offset optimization, Quartic Tensor Power atom iterations, multi-rank block splitting failure analysis, and full 50,880-shot host verifications on column truncations), and the complete script/provenance inventory.
2. [**`chemistry_optimization_conclusions.md`**](chemistry_optimization_conclusions.md): Deep-dive into Benchmark III (Chemistry: DFTHC-SOS). Details the full 50,880-shot host verification results for $C=347$ and $C=345$, the mathematical anatomy of the LCU normalization penalty vs. Toffoli savings, and why greedy atom substitution disrupts joint L-BFGS correlations.
3. [**`CHEMISTRY_ATOM_GENERATION_RESEARCH.md`**](CHEMISTRY_ATOM_GENERATION_RESEARCH.md): Algorithmic design and formulation of Quartic Tensor Power atom generation and Straight-Through Estimator (STE) Givens angle lattice optimization.
4. [**`SUBMISSION_AND_SCALING_ANALYSIS.md`**](SUBMISSION_AND_SCALING_ANALYSIS.md): Mathematical scaling analysis decomposing the geometric mean score $S = (I \cdot D \cdot C)^{1/3}$, proving why chemistry alone cannot reach 14–17×, and outlining the quantitative requirements for 16.6× on DQI.

---

## 2. Benchmark Summary & Current Deliverables

The intact master submission package is packaged in [`submission.zip`](submission.zip) and verified 100% on the official host verifiers:

| Benchmark | Component File / Directory | Scored Cost | Baseline Cost | Improvement Ratio | Host Verifier Shots | Status |
|---|---|---:|---:|---:|---:|:---:|
| **Task 1: Inverse Square Root** | `invsqrt.kmb` | 190,711.96 | 2,076,774.23 | **10.8896×** | 200,540,928 | **PASS** |
| **Task 2: DQI Reed–Solomon** | `dqi.kmb` | 1,612,737,419,463.02 | 24,439,463,583,672.93 | **15.1540×** | 102,400 | **PASS** |
| **Task 3: DFTHC-SOS Chemistry** | `dfthcsos/` (`submission.npz`) | 757,163,489.50 | 3,355,599,588.32 | **4.4318×** | 50,880 | **PASS** |
| **Overall Geometric Mean** | `submission.zip` | — | — | **9.009623×** | Full-scale | **PASS** |

- **Master Archive**: [`submission.zip`](submission.zip) (SHA-256: `6ebc780a1c8027f26080e47818ea5db63dde58ce9ffbf90f46c4e7411b904800`)
- **Scoring & Verification Log**: [`logs/best_submission_score.json`](logs/best_submission_score.json)

---

## 3. Directory Layout

- **`scripts/`**: Reusable Python and CUDA experimental scripts:
  - `scripts/fast_optimize_offsets.py`: Analytical $O(N^2)$ Scipy L-BFGS-B solver for the chemistry SOS gap offset (4.4s).
  - `scripts/test_checkpoint.py`: Automated pipeline running quantization, circuit synthesis, and host verification.
  - `scripts/prepare_cand.py`: Universal candidate preparation script for arbitrary column truncation $C$.
  - `scripts/test_power_atoms_on_incumbent.py`: CUDA Quartic Tensor Power atom generator.
  - `scripts/test_atom_replacement.py`: Atom pruner and replacement evaluator.
  - `scripts/run_atom_replacement_experiment.py`: Atom substitution experiment.
  - `scripts/package_and_score_best.py`: Hermetic packaging and verification script.
- **`runs/`**: Verified candidate outputs and logs:
  - `runs/c347_eval/verify.json`: 50,880-shot host verification of Candidate $C=347$ (PASS, 765.2M).
  - `runs/c345_eval/verify.json`: 50,880-shot host verification of Candidate $C=345$ (PASS, 768.4M).
  - `runs/c347/factor_16k.npz`: Trained factor coordinates for $C=347$.
  - `runs/c345/factor.npz`: Trained factor coordinates for $C=345$.
- **`dfthcsos/`**: Master chemistry submission files (`submission.npz` and `angle_k_*.kmb`).
- **`logs/`**: Host verifier logs and JSON scoring outputs.
