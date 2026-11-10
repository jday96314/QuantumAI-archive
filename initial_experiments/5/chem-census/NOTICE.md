# Attribution

The immutable native DFTHC and kickmix infrastructure is supplied under
Apache-2.0; see `/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/LICENSE.txt`.
No host source, scoring, physical bounds or confidence calculations are changed.

Read-only reused research code:

- `initial_experiments/5/chem-qfit/core.py`: QR PSD/rank-C block projection,
  exact integer-word decoder and canonical research NPZ convention.
- `initial_experiments/5/chem-audit/audit_core.py`: fresh native signed target,
  exact physical quantities, damped strict-gap projection.
- `initial_experiments/5/chem-qfit/greedy_prune.py`: covariance leverage deletion.
- `initial_experiments/5/chem-audit/integer_refine.py`: coupled exact integer
  proposals, adapted locally with output paths redirected inside this directory.
- `initial_experiments/5/chem-qfit/physical.py`: AGY offset derivative and
  integer one-body solver, originally generalized from round-4 AGY.
- `initial_experiments/4/codex/angle_gauge_probe.py`: signed Givens encoder.
- `initial_experiments/5/chem-qrom/{hybrid,build,quadratic,erase,linear}.py`:
  feature row grouping, exact nonlinear GF(2) bank interpolation, streamed
  measurement phase repair, and arbitrary-input XOR lookup compilation.

Imports remain read-only. Every new checkpoint, log, cache, circuit and NPZ is
written inside `initial_experiments/5/chem-census`. Physical/code artifact hashes
are retained in the local manifests; no prior selection is edited.
