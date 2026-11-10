# Round 5 chemistry assumption audit

**Complete archive:** [submission.zip](submission.zip). Full host-metric PASS; local combined score **9.4712493866**. Exact subtree-pruning saving is **0 T**, so chemistry remains banked. DQI uses the newly accepted trimmed decoder.

Assigned directory: `initial_experiments/5/chem-audit`. All outputs stay here;
imports of host and earlier experiments are read-only. Two CPU threads per
process; GPU jobs run sequentially. Full native runs use the shared
`/tmp/r6-fullverify.lock`.

Preserved starting selection: round-4 Codex
`candidates/norm_cycle2_signed`, R1/B719/C350, eight bits, J7, chi10.
Reference chemistry 708,547,388.9553; inverse 182,956.5272; DQI
1,598,204,229,385; combined reference 9.3677221531. **Phase A replacement
native verified: 708,156,073.7876**, selected in local `selection.json`.

**Phase B stopped after the requested frontier audit.** The bound is universal; no tested scratch fit beats the incumbent error/Lambda point. The prepared noise/descent/growth branch remains unrun under this stopping rule.

## Hypotheses and protocol

Phase A: the discrete incumbent improvement frees about 0.48 mHa, which may
permit lower normalization via constrained weight/offset optimization.
Start from the actual selected artifacts and the qfit discrete control;
keep the physical generator gauge aligned and the native fixed shift.
Only a full default native pass with lower actual `lambda_eff*T*Q` can replace
the selection.

Phase B: audit whether the normalization lower bound is universal, then fit
the error/Lambda frontier from scratch on R1/B719/C350, R1/B639/C350,
R2/B447/C175, R2/B479/C175, and an aggressive host-exact frontier shape.
Use the fresh native FCIDUMP and fixed shift, signed full pair metric,
residual + Lambda penalty + gap term (Low et al., arXiv:2502.15882 Eq. 41).
At checkpoints record float, direct rounding and exact integer refit errors
and Lambda. Expected-error angle noise uses sampled uniform half-LSB noise,
without a straight-through estimator; cap norms and Lambda. Also test small
integer-native growth. Stop Phase B if the bound is universal and no tested
regularized fit improves the incumbent error/Lambda point.

## Attribution

Read-only reused code: round-4 Codex `angle_gauge_probe.py`,
`build_selectswap_chem.py`, `monomial_unlookup.py`; round-5 chem-qfit
`core.py`, `physical.py`, `discrete.py`, `oracles.py`; AGY's fast offset
derivative as generalized in `chem-qfit/physical.py`. Host DFTHC math and
circuits are Apache-2.0: `../../../from_host/verifiers/LICENSE.txt`.
No verifier, scoring, confidence, host or prior artifact is modified.

## First numerical checkpoint

`results/calibration.json` independently reproduces native error2
**19.9953044905532 mHa**, Lambda **303.5283147431962**, gap
**3.915157472365536**, lambda_eff **48.59417497035283**. The freshly native
pair-space error differs by less than 1e-15 Ha. The older spectrum's target
differs from the fresh native target by 3.1116404561e-7 Ha Frobenius; all new
optimization uses the fresh target. Source file SHA256 hashes are saved.

The qfit discrete dictionary/weight control gives **19.516374182752 mHa**,
but its weights in the aligned original physical gauge initially have Lambda
**403.99325** and gap **4.99487** after an attempted Newton gap projection.
This is a numerical accuracy improvement, not banked cost margin. A damped
line search replaces the undamped projection so an unreachable target does
not produce oscillating offsets. No success is inferred from projection.

### Bound audit

Located in round-4 Codex `NORMALIZATION_BOUND.md` and
`results/normalization_global_bound.json`. It bounds **every native
factorization at the verifier's fixed shift**, not just the selected model.
It assumes real finite factor weights/offsets, normalized projector rows,
the native LCU one-norm and SOS definitions, fixed S/beta1/ecore/electron
count/upper energy, the strict gap floor, and the physical gap domain. It
does not assume a shape, rank, angle grid, dictionary, factor gauge, target
residual, or closeness to the incumbent.

More precisely, the numerical bound needs the native value of **c**, not a
particular S matrix: S enters only h1_prime, whose positive-eigenvalue sum is
nonnegative for every real symmetric h1_prime. Changing S alone with c held
fixed would not evade the inequality. Changing beta1 (or the other constants
in c), the gap floor, or the native normalization/SOS definitions can change
the numerical bound; those are fixed by this verifier.

Let c = E_upper - ecore + eta*beta1 = -559.1324731485038. For each factor
column t=sum(w), l=sum(abs(w)), and offset W, l>=abs(t). Direct expansion:

```
u = 2 Lambda - gap
  = -c + 2 sum(max(eig(h1_prime),0))
       + sum[(t*t+l*l)/2 + W*t + l*abs(W+t)].
```

Each bracket is convex in W, minimized at W=-t, where its value is
(l*l-t*t)/2>=0. Thus u>=-c and
lambda_eff^2=gap*u > 3.914898*559.1324731485038. The strict universal bound
is **46.78617959252637**, limiting normalization-only improvement to
**1.03864379168x**. The independent native identity differs by 2.5e-12;
the smallest incumbent bracket is positive (0.0016104531). The bound places
no restriction on circuit/resource improvements or the combined score.

At the **fixed, measured** incumbent resource budget T=30,504 and Q=478,
this normalization bound implies cost >682,185,167.45. Therefore a 567M
result requires a lower resource product; a weight/offset refit alone cannot
reach it. This statement holds that resource budget fixed and is not a bound
on all circuits for the shape.

The requested objective is the residual/epsilon_reg + Lambda/lambda_reg +
ReLU((gap-E_reg)/E_reg) expression in Low et al. In the available
[arXiv v1 HTML](https://arxiv.org/html/2502.15882v1#S4.SS1), it is numbered
**42**; equation 41 defines the symmetry-shift family. Our shift is held
fixed to the verifier rather than optimized over that family.

## Reproduction

From this directory: `source env.sh`, then
`"$AUDIT_PY" calibrate.py > results/calibration.log 2>&1`.

## Phase A: banked native improvement

The 8,000-step unconstrained-error seed refit and a 30,000-step continuation
reduce Lambda 404 -> 307 but do not beat the incumbent. Their saved feasible
best has lambda_eff about 48.876. The failure concerns this initialization.

Projecting the **incumbent's physical generators** into the discrete words
instead starts at error2 19.982566 mHa, Lambda 304.29997, weight norm 3.72319.
The next 30,000-step exact-native constrained refit gives error2
**19.9899725541 mHa**, Lambda **303.2039440**, gap **3.9150496064**,
lambda_eff **48.5673692**. A gauge-only 30-sweep Jacobi rotation preserves
the augmented coefficient Gram to 3.4e-13 and reduces Lambda to
**303.2035505859**, lambda_eff **48.5673374743**.

The analytic offset SLSQP solver succeeds in 259 iterations but its proposed
gap 3.91515 raises the product slightly; it is rejected. At differentiable
probe points its gradient directional errors are <2.5e-7 for gap and <1.5e-4
for the squared product. Direct finite differences at near-zero |w_B| kinks
were initially rejected; checks now use explicitly perturbed differentiable
points. Signed re-encoding is required after integer moves; it preserves
every projector to 1e-12 before circuit construction. A temporary helper
named `select.py` shadowed Python's standard module and was renamed
`audit_selection.py` before any successful build.

One-body descent reaches **11.99537293 mHa**. All eleven compiled ANGLE
oracles obey the 128-qubit cap. Static host resources remain T **30,504**,
Q **478** (141 persistent + 337 temporary). Full default native verification
under the shared lock gives **PASS**, **50,880 shots**, actual score
**708,156,073.7876049**, agreeing with the static estimate. Both errors and
all strict floors pass. Canonical digest:
`cafc58c606bdc7393f08b74b0a50789b0d5219edf502ebf352bc09d7da343c00`.

`audit_selection.py` independently binds the circuit/NPZ bytes to the native
digest, checks the default shot count, physical accuracy and actual objective,
then writes `selection.json` with all twelve artifact SHA256s. Chemistry
improvement **1.00055258x** (about **0.0552%**); the 567M and 420M targets
remain unmet. The original round-4 candidate and its combined reference are
untouched; inverse and DQI are unchanged. Local replacement directory:
`candidates/banked`.

A second full-default native run with the predeclared independent salt
`round5-chem-audit-bank-independent` also passes all **50,880 shots** with
the same digest and score. This checks the measurement-based erasure on an
independent simulator stream; no salt/hash search is performed.

## Phase B implementation and finite search

`frontier_fit.py` starts with random unit vectors and Gaussian factor weights
and trains Low et al.'s norm-residual + Lambda + ReLU gap objective. S and
beta1 stay fixed; offsets are optimized via the primitive identity weights.
Native exact physical quantities use float64. The first scratch run uses
float64 throughout. Subsequent fast runs may compute the **full, direct**
residual in float32 (no truncated spectrum or PSD clipping), then recompute
float/rounded/refit checkpoints in float64 against the fresh native target.
Their training-vs-checkpoint squared-error discrepancy is explicitly logged;
all potential acceptances still require full native verification.

Requested shape coverage: R1/B719/C350, R1/B639/C350, R2/B447/C175,
R2/B463/C175, R2/B479/C175, R4/B180/C88. The last is the aggressive frontier
shape: host-exact PREPARE/Givens plus realizable eight-bit ANGLE upper bound
T=18,391, Q=443, conditional cost **395,907,094.04** at the original
lambda_eff. This is a resource estimate without a physically passing fit.
`resources.py` also covers seven-bit choices and records conditional costs
at the universal normalization bound; these are not universal cost lower
bounds, because the ANGLE count is an upper bound.

Each shape gets a scratch seed at lambda_reg=300 followed by explicitly
labeled continuations at lambda_reg=30 and 3000. Thus the regularization
sweep shares a scratch-trained basin; it is not presented as independent
random restarts. Initial R1 controls cap weight Frobenius norm at 5; the
broader multiblock and regularization sweep cap it at 100 so this does not
exclude any Lambda<=320 solution (the native Lambda formula already implies
norm(w)<=2*sqrt(Lambda)<36). Robustness runs use the stricter cap of 5. Exact
unconstrained diagnostic refits may exceed this and have their norms and
Lambda recorded, rather than being silently treated as bounded candidates.

`integer_refine.py` reuses qfit's coupled integer proposals but loads the
fresh native metric. The bounded mode accepts exact refit proposals only
when weight norm<=5 and Lambda<=320; an explicit unbounded control can
measure the normalization tradeoff. Offsets and all three checkpoint states
are saved. No straight-through estimator is used.

The initial float32 B639 preflight stopped at a 2.82e-7 squared-error
discrepancy. That failed assertion is preserved in
`results/failed_precision_guard_b639.*`. The guard now measures the error
**norm** difference directly and requires <1e-5 Ha (0.01 mHa); exact native
float64 checkpoints and acceptance thresholds are unchanged. This is a
numerical preflight correction, not evidence against the physical design.

Plotting uses Matplotlib installed only into local `plot_deps/`; the native
verifier/training environment is unchanged. `plot_frontier.py` is a separate
plot-only process reading saved JSON, and exports PNG and SVG.

The first B463 GPU attempt ran out of shared device memory before completing
its first checkpoint; `results/failed_gpu_memory_b463.log` preserves it. A
20-step CPU preflight completed and confirmed the algebraic one-body
contraction at 2.3e-13. CPU training is too slow for the main sweep. Subsequent
training keeps exact float64 checkpoint/refit work on CPU to lower GPU memory
use. The restarted B463 GPU run progresses independently from the failed
attempt; its missing-shape fit is not counted until it finishes.

Operational scope correction: the initial AGENTS.md filename scan filtered
excluded paths instead of pruning traversal. It returned no matches, and no
excluded file contents were opened. Subsequent repository reads used explicit
allowed paths. Phase A/B results use no data from the excluded concurrent
directories. The later packaging follow-up reads the DQI acceptance evidence
and selected artifact under the user's explicit authorization for that step.

Commands (all from this directory after `source env.sh`):

```bash
"$AUDIT_PY" refit.py --fit parameters/discrete_seed.npz --out parameters/banked_margin.npz --steps 8000 --lr 2e-6
"$AUDIT_PY" refit.py --fit parameters/banked_margin.npz --out parameters/banked_margin_long.npz --steps 30000 --lr 2e-6
"$AUDIT_PY" transfer_seed.py
"$AUDIT_PY" refit.py --fit parameters/transferred_seed.npz --out parameters/transferred_refit.npz --steps 30000 --lr 2e-6
"$AUDIT_PY" gauge.py --fit parameters/transferred_refit.npz --out parameters/banked_gauge.npz --sweeps 30
"$AUDIT_PY" build_candidate.py --fit parameters/banked_gauge.npz --tag banked --offset-steps 400 --build
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 flock /tmp/r6-fullverify.lock "$AUDIT_PY" "$AUDIT_ROOT/from_host/verifiers/dfthcsos/verify.py" candidates/banked --salt round5-chem-audit-bank
"$AUDIT_PY" audit_selection.py --candidate candidates/banked --result results/banked_native.json
```


<!-- baseline-start -->
## Measured scratch baseline checkpoints

Latest/terminal checkpoints at lambda_reg=300; errors are mHa. These are
numerical two-body fits, not native passes. Other physical gates still apply.

| Shape | Step | Float error | Rounded error | Refit error | Float Lambda | Float gap | Refit Lambda |
|---|---:|---:|---:|---:|---:|---:|---:|
| R1/B719/C350 | 20000 | 30.2700 | 203.1514 | 137.6643 | 471.7476 | 3.863399 | 572.4446 |
| R1/B639/C350 | 20000 | 29.8641 | 211.0735 | 146.5297 | 466.9635 | 3.893561 | 535.1368 |
| R2/B447/C175 | 20000 | 31.3639 | 193.7871 | 132.0498 | 432.5435 | 3.867388 | 508.8856 |
| R2/B463/C175 | 20000 | 30.3886 | 190.2596 | 131.4224 | 420.4702 | 3.851475 | 503.0335 |
| R2/B479/C175 | 20000 | 29.4405 | 175.5948 | 124.6285 | 453.4937 | 3.811888 | 523.5530 |
| R4/B180/C88 | 20000 | 46.9970 | 202.0299 | 155.5963 | 408.6004 | 3.845036 | 428.1509 |

Reproduction: the first R1/B719 command is recorded below; all remaining
baseline and regularization jobs are in `run_scratch.sh`, which skips
completed baseline runs. `results/code/` preserves trainer/helper versions
for the later runs. JSON arguments and NPZ hashes identify measured states.

```bash
AUDIT_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$AUDIT_PY" frontier_fit.py --R 1 --B 719 --C 350 --seed 11 --lambda-reg 300 --weight-cap 5 --steps 20000 --checkpoint 2000 --precision float64 --checkpoint-device same --tag scratch_r1b719_l300
bash run_scratch.sh
# Prepared but unrun under the Phase B stop rule: run_robust.sh
"$AUDIT_PY" summarize_frontier.py
"$AUDIT_PY" screen_physical.py
```

<!-- baseline-end -->


<!-- final-findings-start -->
## Final Phase B findings and stopping decision

**Stopped under the user-defined condition.** The 46.7861796 bound is
universal for the native fixed-constant problem, and none of the tested
from-scratch Lambda-regularized float fits reaches error2 <=20 mHa with
Lambda below the incumbent 303.5283147. There are six random scratch
starts (20,000 steps each) and twelve regularization continuations
(6,000 steps each). The continuations share those six basins and follow
20-sweep column-gauge reductions; they are not independent random starts.

| Shape | Best float error (mHa) | Lambda there | Best float error at Lambda <=320 | Lambda there | Rounded error there | Refit error there | Refit Lambda |
|---|---:|---:|---:|---:|---:|---:|---:|
| R1/B719/C350 | 28.3892 | 388.5595 | 32.9517 | 305.3164 | 202.3643 | 137.4017 | 406.8031 |
| R1/B639/C350 | 27.8265 | 375.1661 | 32.3612 | 305.3089 | 198.6241 | 141.1662 | 381.5782 |
| R2/B447/C175 | 29.7133 | 347.6691 | 33.8375 | 300.3409 | 184.3911 | 128.8395 | 381.3838 |
| R2/B463/C175 | 28.8388 | 349.0267 | 33.2078 | 301.3119 | 189.8069 | 129.6673 | 396.9761 |
| R2/B479/C175 | 27.9644 | 360.0598 | 32.4306 | 302.3705 | 185.5818 | 129.0490 | 377.1894 |
| R4/B180/C88 | 45.8619 | 327.9249 | 48.7567 | 298.6579 | 198.8383 | 153.8508 | 320.0866 |

The scratch-family minimum rounded/refit error over all saved checkpoints is **123.6072 mHa**.
The complete CSV has **114 checkpoint triples / 342 metric rows**.
Every triple records float, direct rounding and exact rank-C refit error,
Lambda, gap, effective normalization and weight norms. Saved research
NPZs, JSON arguments and source snapshots identify the exact measured states.

The maximum recorded float32-vs-float64 error-norm discrepancy is **0.001888 mHa**;
the maximum recorded physical-contraction discrepancy is **3.41e-13 Ha**.
These are numerical consistency checks, not native acceptance. The final
`physical_shortlist.json` is empty: no scratch integer state meets the
two-body accuracy gate, even before building one-body words and oracles.

The expected-half-LSB-noise, integer-descent and integer-growth branch was
**prepared but not executed**, because the explicit Phase B stop condition
was met after the frontier audit. `run_robust.sh`, `integer_refine.py` and
`integer_grow.py` are unrun research implementations, not evidence of
success or failure for those methods. The CPU preflight is an implementation
check only and is not counted as a completed frontier fit.

The negative fit results apply to these seeds, budgets, parameterization
and regularization paths. They do not establish global infeasibility of
these shapes or a universal chemistry-cost bound. Indeed the preserved
known integer basin is feasible at R1/B719; the random searches did not
reproduce it. Other basins, gauges, longer optimization and oracle designs
remain possible. The rigorous universal result concerns normalization only.

**Final selection remains `candidates/banked`, cost 708,156,073.7876, improvement 1.000552583x.**
It passes two full-default native salts with 50,880 verified shots each.
The <=567M target and <=420M stretch were not reached. No later research
fit replaces the native-verified margin gain.

Artifacts: [all checkpoint values](results/frontier_checkpoints.csv),
[frontier plot](results/frontier.png), [vector plot](results/frontier.svg),
[native result](results/banked_native.json), [independent native result](results/banked_native_independent.json),
[selection](selection.json), [source/artifact hashes](artifact_manifest.json).

<!-- final-findings-end -->


<!-- package-followup-start -->
## Follow-up: exact subtree probe and complete package

**Probe outcome: zero projected Toffoli saving, including optimized address order.**
Every one of the 795 full seven-word XOR-difference payloads is distinct
in every oracle. Every one of the 315,615 distinct-input pair XORs is also
distinct in every oracle. Thus a row permutation cannot create a constant
subtree containing two valid payloads, nor repeated pair-difference leaves.
The selected emitter already skips empty ranges and collapses constant
subtrees; the final artificial padding row is already handled.

The exact selected-emitter count is 396 lookup + 49 low-bit selection +
37 ANF erasure = **482 Toffolis per oracle per pass**. The fixed 38
one-body pairs activate all 49 bank bits, so permuting dictionary rows
does not reduce that component either. The counter equals the loaded
accepted KMB operation statistics for all eleven oracles.

| Oracle | Current T/pass | Current round-trip T | Current-order saving | Optimized-order saving |
|---|---:|---:|---:|---:|
| 0 | 482 | 964 | 0 | 0 |
| 1 | 482 | 964 | 0 | 0 |
| 2 | 482 | 964 | 0 | 0 |
| 3 | 482 | 964 | 0 | 0 |
| 4 | 482 | 964 | 0 | 0 |
| 5 | 482 | 964 | 0 | 0 |
| 6 | 482 | 964 | 0 | 0 |
| 7 | 482 | 964 | 0 | 0 |
| 8 | 482 | 964 | 0 | 0 |
| 9 | 482 | 964 | 0 | 0 |
| 10 | 482 | 964 | 0 | 0 |
| **Total** | **5302** | **10604** | **0** | **0** |

Savings columns include both forward and reverse passes. All fourteen
tested common permutations give the same count: current order, eleven
per-oracle lexical orders, a combined lexical order, and a common
nearest-neighbour/2-opt order. One-body rows stay fixed; dictionary U,
angle-word and weight rows travel together, with exact inverse-permutation
array checks. The equality-collapse objective is identically zero, so the
zero optimum is certified independently of the heuristic search.

The secondary Hamming path objective improves from 193,950 to 169,293
(12.7131%); this does not save Toffolis in this emitter.
Chemistry T remains **30,504**. The 3% build threshold is 915.12 T, hence
at least **916 integer Toffolis**. Projected saving is **0 T / 0%**, so
**no pruned family was built or native-verified**. This certificate is
specific to exact whole-payload equality under permutations. It does not
bound bitwise factoring, alternate gauges, linear transformations or
other oracle architectures.

### Packaged accepted artifacts

- Chemistry: unchanged `candidates/banked`, canonical digest
  `cafc58c606bdc7393f08b74b0a50789b0d5219edf502ebf352bc09d7da343c00`.
- Inverse: round-4 Codex `candidates/zero_head_affine.kmb`.
- DQI: `initial_experiments/5/dqi-core/candidates/trimmed_decoder.kmb`.
  The acceptance marker appeared while packaging was being prepared.
  The first assembly guard stopped before writing an archive; the updated
  default/independent/paired full native reports were then checked.
  All have 102,400 shots. Paired runs have identical failure counts, R
  and Q, with exactly 25,596 fewer mean Toffolis. The default accepted
  report is the reference; the lower independent-salt cost was not selected.

Archive: **[submission.zip](submission.zip)**, **110,392,317 bytes**, exactly 14 members.
SHA-256: `3b5805520a850216f0b7a42fa2742007777f2c841e27df3dcdea1ddb2a50810a`.

**Full unmodified host-metric validation PASS.** `from_host/metric.py`
read/validated every circuit and required NPZ array, ran all three native
verifiers with full defaults and verification_shot_divisor=1, then computed
`metric.score`. CPU affinity was restricted to two CPUs; native worker
defaults therefore used at most two workers. Only the local input-directory
configuration was set; no verifier, confidence or scoring function changed.

| Task | Verified shots | Measured cost in full metric run |
|---|---:|---:|
| invsqrt | 200,540,928 | 182956.466084 |
| dqi | 102,400 | 1.54722161865e+12 |
| dfthcsos | 50,880 | 708156073.788 |

**Measured local combined score: 9.4712493866031.**
Accepted-report reference combined score: **9.4455168217812**.
These are deliberately separate: DQI sampling failures and the host
Clopper–Pearson repetition factor vary with the salt. The metric run used
one predeclared fixed local salt, with no seed/hash/salt search:
`chem-audit-complete-archive-validation-2026-10-24-fixed-v1`.

The archive was independently reopened after verification. Every member
was read to EOF for CRC/content SHA checks; inverse/DQI circuit hashes,
chemistry canonical digest, full shot counts, task costs and combined
score all match the fresh native reports. The original accepted artifacts
and all round-4 inputs remain untouched.

Packaging layout/hash binding adapts read-only round-4 Codex
`package_submission.py` and `audit_submission.py`; the metric and native
Apache-2.0 source stay unchanged under `from_host/`. DQI attribution and
its accepted full-run evidence are retained in
`results/dqi_core_selection_evidence.txt` and `package_trimmed_dqi_proofs.json`.

```bash
source env.sh
"$AUDIT_PY" subtree_probe.py > results/subtree_probe.log 2>&1
"$AUDIT_PY" assemble_archive.py > results/package_assembly.log 2>&1
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 flock /tmp/r6-fullverify.lock "$AUDIT_PY" package_metric.py > results/package_metric.log 2>&1
"$AUDIT_PY" package_audit.py > results/package_independent_audit.log 2>&1
"$AUDIT_PY" record_package_results.py
```

Evidence: [probe JSON](results/subtree_probe.json), [all permutation/oracle counts](results/subtree_probe.csv),
[assembly](results/package_assembly.json), [fresh metric result](results/package_metric.json),
[independent ZIP audit](results/package_independent_audit.json), [complete summary](submission.summary.json).

<!-- package-followup-end -->
