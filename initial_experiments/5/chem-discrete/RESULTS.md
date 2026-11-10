# Round 5 continued: direct discrete chemistry

**Final selection: no replacement.** The 1.10x, 1.20x and 1.35x targets were
not reached. The best new native-verified model costs **727840914.114845**,
which is **2.722969% worse** than the selected incumbent. Improved integer
dictionaries, capacity certificates and native JSONs are preserved.

The selected incumbent costs **708547388.9553173**. Its full native reproduction passes with
50880 shots, T=30504, Q=478, error1=11.9863154428 mHa,
error2=19.9953044906 mHa, Lambda=303.5283147432,
E_gap=3.9151574724, lambda_eff=48.5941749704.
Evidence: `results/incumbent_native.json` and stderr.

All writes are confined to this directory. The completed chem-qfit experiment,
earlier selected artifacts, and native host implementation are read-only.
`core.py`, `physical.py`, and `oracles.py` start from the attributed chem-qfit
implementations; imported signed encoding and binary column-gauge optimizer
are read-only round-4 Codex utilities. New search code does not train
continuous orbital vectors or round continuous fitted dictionaries.

## Controls and search methods

`reproduce.py` binds source checkpoint hashes, decodes actual words, measures
the full signed symmetric-pair tensor residual, and preserves 0/1/4/20/100
sweep controls. Final native-target checks reproduce the saved split at
22.6232801615 mHa and the 100-sweep state at 22.5990813893 mHa. The accepted
B719 dictionary reproduces exactly; its fixed-word refit is 19.5330305538 mHa
and the saved changed-word control is 19.5163741828 mHa.

* `joint_refit.py`: matched fixed-word PSD rank-C control, using a QR basis
  containing both block spans and joint whitened factor L-BFGS. Its exact
  compact loss includes an explicitly computed outside-span residual and is
  checked against the full tensor. It diagnoses block solver basins.
* `search.py`: finite integer multi-angle row moves, signed residual Ritz-mode
  directions, full cross-block finite-feature interactions, and binary
  annealing for simultaneous multi-vector proposals. Every accepted proposal
  has a full signed residual after separate PSD rank-C block refits. Thermal
  mode also accepts actual uphill integer states and restores the best state
  between temperature cycles. Single mode supplies a matched simpler control.
* `replace.py`: entirely new representable residual-mode atoms, quantized
  dictionary combinations, exact-word donor atoms, multiple replacements,
  and correlated donor/receiver swaps between blocks. The common-span proxy
  only orders trials; the full PSD rank-C model decides acceptance. Actual
  uphill states are allowed at a finite integer temperature.
* `normalize.py`: multi-block weights and identity offsets, native fixed
  shifts and actual SOS gap, target floor margin 0.000402 Ha and accepted
  strict-floor margin at least 0.000152 Ha.
  An augmented stage approaches feasibility while optimizing normalization;
  the final `--feasible` stage uses backtracking and accepts only feasible
  one-body/two-body/gap/Lambda iterates. It preserves best feasible states.

Every search keeps its best checkpoint and a JSON trajectory beside the NPZ;
periodic and improved snapshots live in `trajectories/`. Proposals do not
discard superior intermediates when exploring uphill moves.

## Stronger fixed-dictionary capacity results

`capacity.py` solves the arbitrary-symmetric-coefficient block relaxation using
principal angles between the two block spans. It excludes cross-block tensor
coefficients absent from the submitted model. This is stronger than the union
span relaxation. Final analyses preserve immutable `capacity_snapshot_*.npz`
files and hashes; the earlier exploratory B639 calculation used a moving
checkpoint and is superseded by `results/capacity.json`.

* Original saved split B479: joint-span relaxation approximately **17.783814 mHa**;
  block-linear relaxation **19.4857969641 mHa**.
* Saved changed-word Claude B479 control: joint-span approximately **17.849205 mHa**;
  block-linear **20.3017634201 mHa**. These particular words cannot pass,
  even without PSD or rank restrictions.
* Saved B447 control: joint-span approximately **20.596204 mHa**;
  block-linear **23.3046923853 mHa**. Dictionary changes are necessary.

More decisively, the original saved split B479 has a tight **PSD** capacity
certificate, allowing unlimited rank inside each block:

```
native dual lower bound = 20.156389669428176 mHa
native primal residual = 20.156416404145453 mHa
max eigenvalue of each dual block compression < -9.99e-12
```

Thus those exact words cannot reach 20 mHa for any C. This is a fixed-dictionary
obstruction, not a numerical local-minimum claim and not a shape-wide theorem.
`admm_refit.py` solves the convex unlimited-rank PSD problem; the nonconvex
rank-175 ADMM variant was unstable and did not improve its starting candidate.
`psd_certificate.py` constructs a feasible dual residual D. For every model
matrix M in the sum of the two PSD block cones,
`||G-M||^2 >= 2<G,D> - ||D||^2`, because each `Q_r.T D Q_r` is negative
semidefinite. The final scalar safety shift is explicit. The full dual matrix,
words and JSON are retained in `candidates/admm_fullPSD_B479_polish_dual.npz`
and `results/admm_fullPSD_B479_polish_dual.json`. This is a float64 numerical
certificate with substantial margin over the calibrated native/pair discrepancy.

The subsequent `psd_geometry_B479_full` search changes integer words using
the convex PSD residual. It has already beaten the original dictionary's
certified floor. Every twentieth checkpoint receives a separate actual
C175 PSD/rank refit; these audits remain above threshold. An improved relaxed
geometry is not reported as a passing C175 factorization.

## Native passing normalization controls

| Candidate | error1 mHa | error2 mHa | Lambda | E_gap | lambda_eff | T | Q | Native cost |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Selected incumbent | 11.986315 | 19.995304 | 303.528315 | 3.915157 | 48.594175 | 30504 | 478 | 708547388.9553173 |
| Feasible B711 | 14.588863 | 19.991163 | 322.280780 | 3.915515 | 50.084522 | 30408 | 478 | 727979729.6680878 |
| Gauge-escaped B711 | 14.588863 | 19.991163 | 322.223489 | 3.915515 | 50.080043 | 30408 | 478 | 727914625.9058473 |
| Final feasible B711 | 14.477221 | 19.991204 | 322.122686 | 3.915958 | 50.074972 | 30408 | 478 | 727840914.114845 |

All three new B711 directories pass full 50368-shot native forward/reverse
verification, strict floors, both residual bounds and the 128Q oracle cap.
Their JSONs and exact submission directories are preserved. None beats
the incumbent; all three are explicitly unselected. The final stage uses
exact native eigenvector matching and retains feasibility throughout.

The final feasible optimizer checks the native Hungarian one-body matching
at every trial. An earlier control used a fixed matching; its reported native
metrics above are authoritative. Configured target gap is 3.9153, and accepted
final iterates require at least 3.91505, a **0.000152 Ha** strict-floor margin.
Both native physical bounds remain satisfied throughout the final stage.

The new physical gradient agrees with autograd to **2.82e-18** and finite
differences to **5.04e-14**. Rank-aware exchange spectral screens agree with
direct new-dictionary eigensolves to **1.43e-14**. See their control JSONs.
All-address arbitrary incoming-target XOR/reverse tests pass for the new
B711 and B479 oracle families, including the 127Q padded B479 construction.

## Resource and normalization results

Actual signed-oracle resource measurements give B479/C175 T=27981, Q=461,
break-even lambda_eff=54.9293861; B447/C175 T=26171, Q=461,
break-even=58.7283310; B639/C350 T=27794, Q=478,
break-even=53.3322556. The required multipliers for a 1.10x win are these
numbers divided by 1.10. Full C/shape sweeps are in `results/rank_frontier.json`.
Increasing C on B479 reaches about 20.73 mHa under short cyclic refits but
rapidly worsens T/Q; converged full PSD refitting supplies the stronger bound
above. Neither experiment establishes a general B479 shape floor.

The B639 integer dictionary improved from **102.0062 to 34.0700 mHa** in
160 steps, while its actual lambda_eff is **103.7044** and conditional resource
product **1377768456.83**. The normalization-aware continuations finish at
**29.8437484886 mHa** after **802 total integer iterations** (160 + 444 + 198).
Normalization was measured during the search, before any physical pass.
Final raw lambda_eff is **103.5701318**; twelve column-gauge sweeps reduce it
to **80.8870924**, preserving the physical Gram matrix within 2.71e-13.
Actual T/Q remain 27794/478; the conditional product is
**1074628054.2268462**, still uncompetitive and physically invalid.
This is a finite normalization screen, not a normalization lower bound.
Whole-atom, donor/receiver, subtree, and true
uphill searches retain their best initial state when no improvement survives.

## Early normalization finding

The saved B479 split has native error2=22.6232801615 mHa and
lambda_eff=74.8978991713 after four column-gauge sweeps. The B671 integer
control has error2=21.2258890529 mHa and lambda_eff=93.0436587963 under the
same screen. Their gaps are about 4.08 Ha; both fail accuracy and cost.
These are finite gauge/offset screens, not normalization lower bounds.
Complete records: `results/start_legacy_exact_split_B479_physical.json` and
`results/start_anneal_B671_physical.json`.

## Matched direct-native comparison

All four controls start from the identical `refit100_legacy_exact_split_B479`
checkpoint. Searches use seed 101 and the direct native-reference pair target.
The requested 650-eigensolve limit rounds to **652 actual block eigensolves**,
because a trial executes complete two-block, three-cycle refits. The fixed
control uses exactly 326 cycles = 652 block eigensolves.

| Method | Final error2 mHa | Integer row moves | Uphill accepts | Block eigensolves | Loop seconds |
|---|---:|---:|---:|---:|---:|
| Fixed words, refit only | 22.5468939677 | 0 | 0 | 652 | 164.40 |
| Single-angle control | 22.5444806203 | 102 | 0 | 652 | 180.14 |
| Cross-block coordinated QUBO | 22.5470165230 | 336 | 0 | 652 | 158.49 |
| Annealing + subtree proposals | 22.5473270161 | 510 | 4 | 652 | 184.82 |

Coordinated and thermal moves **do not beat the simpler control at this
matched spectral-work budget**. Most progress is weight refitting. Timings
are descriptive on a shared GPU; eigensolves match expensive refit work while
proposal overhead differs. Original 1000-second explorations and this serial
native-target repeat are retained. Commands, return codes, hashes and wall
times are in `results/native_matched_manifest.json`; full trajectories are
`candidates/native_matched_*.npz.json`.

## Additional discrete searches and failure classification

* **Best new actual C175 dictionary:** normalization-aware B479 reaches
  native error2 **22.4021236984 mHa**, error1 **14.9389707274 mHa**,
  lambda_eff **73.4511467404**, T/Q **27981/461**, conditional product
  **947464043.5313091**. It fails accuracy and cost. Exact physical state:
  `candidates/evaluated_norm_coupled_B479_merit.npz`.
* **PSD-residual geometry:** 148 integer iterations and 1637 accepted row
  moves lower the old PSD floor to native error2 **20.0671507364 mHa**.
  This still fails and uses C479: T/Q **46984/466**, lambda_eff **104.7832075**.
  Actual C175 audits at every twentieth checkpoint finish at
  **22.5860003988 mHa**. Neither the relaxed residual nor the larger C is
  presented as a passing C175 model. All audited dictionaries are retained.
* **Whole-atom basin changes:** B479 ran 56 iterations/1001 seconds,
  accepted 39 actual uphill moves and returned no superior best dictionary.
  B447 ran 40 iterations/1214 seconds, with 23 uphill accepts and integer
  quenches, also retaining its initial best. Donor/receiver swaps, multiple
  replacements, new residual-mode atoms and quantized combinations were
  screened by rank-C spectral gains and checked after full tensor refits.
  These finite failures do not establish impossibility of different words.
* **Independent initialization:** alternating integer row-support selection
  and continuous factor-column gauges plateaued near **46.3811 mHa** for
  B479 and **235.1558 mHa** for random balanced B447. Actual tensor fits are
  worse than the saved integer starts; the factor-span surrogate was not a
  useful physical initialization predictor.
* **Numerical block basin:** joint whitened rank-175 L-BFGS improves the
  unchanged split dictionary from 22.6233 to **22.4178 mHa** in 300 outer
  steps. Nonconvex rank-175 ADMM is unstable and never beats its input.
  Convex full-PSD ADMM converges enough to supply the dual obstruction above.
  Rank limitations, coefficient-span limitations and numerical local minima
  are explicitly distinguished.
* **B639:** 102.0062 -> 34.0700 -> 30.2548 -> **29.8437 mHa**, through
  exact-word moves and continuous weight/offset refits. This is useful
  discrete progress, not a physical pass or a score improvement. B511 was
  not pursued because easier tested shapes never became physically valid
  and competitive together.

No result establishes that all integer B447, B479 or B639 dictionaries are
impossible. Proved obstructions apply only to preserved exact word tables.
The normalized changed-word B479 basin has a block-linear relaxation below
20 mHa and remains unresolved at actual C175.

## Exact native metric audit

The inherited signed-spectrum target differs from the direct native pair
target by **3.1116404561e-7 Ha** in Frobenius norm. Historical traces used
that inherited target; measured endpoint discrepancies are below one
nanohartree. `native_target.py` now constructs the target directly from the
checksum-validated FCIDUMP, native float32 TEI and fixed float32 shifts.
`core.py` and new searches use this direct target.

The matched comparison was repeated on it, the last B639 continuation uses
it, and the final checkpoint audit recomputes saved numerical states with it.
The incumbent pair error reproduces native error2 within **9e-17 Ha**.
Native PASS scores always came from the unchanged native verifier. The dual
certificate is independently checked against the direct native target in
`results/native_pair_control.json`. Both native- and spectral-residual dual
proposals are tested for cone feasibility and priced against the **native**
target; the stronger feasible multiplier is preserved.

The final feasible stage's one-microhartree error2 acceptance margin also
exceeds the 0.311164-microhartree target discrepancy, so the earlier target
did not invalidate retained feasibility. One-body matching, fixed shifts,
actual SOS gap and full native final verification were checked separately.

Two exploratory batches exhausted GPU memory; their checkpoints and failure
logs are retained. Small rank-screen eigensolves and one-body searches were
moved to CPU, concurrency was reduced, and final comparisons completed.
Identical-projector donor proposals were excluded after they obscured actual
basin changes. Best states survive temporary uphill moves.

## Saved artifacts and verification

* `selection.json`: no replacement and exact native score comparison.
* `results/native_reproduction.json`: native-target reproduction of the
  strongest completed qfit controls, with source SHA256s.
* `results/numerical_summary.json`: direct-native exact-word metrics,
  decoded-vector checks and hashes, including saved trajectories.
  The final audit covers **1859 checkpoints**; maximum decoded-vector
  disagreement up to sign is **3.34e-16**. Fourteen fully assessed saved
  states pass physics; no priced new feasible state beats the incumbent.
* `results/physical_summary.json`: complete significant-checkpoint errors,
  Lambda, gap, effective normalization, host resource screens and status.
* `results/*_native.json`: four full native PASS JSONs; three are unselected
  new B711 controls and one reproduces the selected incumbent.
* `results/native_digest_audit.json`: each native digest matches delivered
  circuit/NPZ bytes. The incumbent digest remains
  `ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7`.
* `results/oracle_forward_reverse.json`: arbitrary incoming-target XOR tests
  on every address in both directions, including 127Q B479 oracles.
* `candidates/` and `trajectories/`: improved exact-word dictionaries, warm
  controls, best states, periodic states, normalization states and dual arrays.
* `artifact_manifest.json`: preserved research artifact hashes.
* [Progress figure](results/progress.png), [exportable PDF](results/progress.pdf).

Invalid models' finite products are **conditional resource products**, not
accepted native task scores. Only native PASS rows carry accepted task costs;
all new ones are worse than the incumbent.

## Reproduction commands

Run here with `source env.sh` (Python, CUDA, and native C++ runtime paths are
specified there). No external service or changed host verifier is used.

To reproduce measurements of the delivered exact artifacts without rerunning
search, run `bash verify_saved.sh`. `bash reproduce_all.sh` contains the
budgeted numerical experiment commands, including matched comparisons,
capacity certificates, continuations, normalization and native verification.
Fresh duration-limited searches can execute different iteration counts with
different GPU contention; delivered checkpoints and native digests bind the
reported outcomes. The fixed-eigensolve comparison is independently repeatable.
Plotting dependencies are local to `results/plot_dependencies`, with versions
in `results/plot_requirements.txt`; no shared Python environment was modified.

```bash
"$DISCRETE_PY" reproduce.py > results/reproduction.log 2>&1
"$DISCRETE_PY" "$DISCRETE_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  "$DISCRETE_ROOT/initial_experiments/4/codex/candidates/norm_cycle2_signed" \
  --salt round5-chem-discrete-control > results/incumbent_native.json
"$DISCRETE_PY" joint_refit.py --init candidates/start_legacy_exact_split_B479.npz \
  --out candidates/joint_B479.npz --steps 300 > results/joint_B479.log 2>&1
"$DISCRETE_PY" search.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/coupled_B479.npz --mode coupled --steps 200 --seconds 1000 \
  > results/coupled_B479.log 2>&1
"$DISCRETE_PY" replace_gpu.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/rankreplace_real_B479.npz --steps 100 --seconds 1000 \
  --trials 4 --temperature .00001 > results/rankreplace_real_B479.log 2>&1
"$DISCRETE_PY" rank_frontier.py > results/rank_frontier.log 2>&1
"$DISCRETE_PY" normalize.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/norm_B479.npz --steps 4000 --seconds 1500 --lr .00003 \
  > results/norm_B479.log 2>&1
```

No winning submission was produced. The three new B711 submission directories
pass both physical thresholds, strict floors, exact tree and XOR words, the
128Q oracle cap, and full default forward/reverse native verification, but
lose on cost. No numerical or relaxed checkpoint is promoted.

Current improvement targets correspond to costs **644133989.9594** (1.10x),
**590456157.4628** (1.20x), and **524849917.7447** (1.35x). None was reached.
All numerical experiment processes are terminal. Source and host hashes are
in `results/reference_provenance.json`; final byte binding is recorded in
`results/native_digest_audit.json` and `artifact_manifest.json`.
