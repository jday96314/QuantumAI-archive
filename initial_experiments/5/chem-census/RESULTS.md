# Chemistry census: routes beyond the incumbent basin

Status: complete. Local verified compiler improvement selected: 696,942,846.9115;
1.25x target not reached. External selections remain untouched.
Baseline: banked chem-audit cost 708,156,073.7876, lambda_eff 48.5673375,
T 30,504, Q 478. Target: 566,524,859.0301 (1.25x).
Only this directory is writable. Prior results and host code are read-only.
GPU was idle and chem-audit's report/process check showed its work terminal.
Numerical jobs use one GPU and at most two CPU threads in aggregate; full
verifiers run under `flock /tmp/r6-fullverify.lock`.

## Phase 1: distinct routes and estimates

The following are conditional resource estimates, not verified submissions.
Exact host threshold mapping and compiled ANGLE counts will update them.

| Rank | Route | Objective estimate | Biggest uncertainty and decisive probe |
|---:|---|---|---|
| 1 | Integer continuation: remove low-leverage rows from the actual banked dictionary, refit in its physical column gauge, repair exact words at each 16-row boundary | B703 about 695M at fixed normalization; B511 about 523.5M, or 1.35x. Each threshold requires lambda increase below roughly 1.8%, until address-width drops help further | Can deleted directions be replaced without leaving the low-Lambda basin? Compare native signed tensor optimum, exact-word repair, and native Lambda at B703 and deeper thresholds. |
| 2 | Integer-native growth: start from a smaller *integer* subdictionary; quantize residual-directed atoms before insertion, retain incumbent physical generators as a warm reference, and stop at a resource frontier | B511/C350 about 523.5M; B575 about 600M, requiring oracle savings or lower C. Goal lambda at B511 is about 52.6 | Does growth recover useful missing span at bounded normalization, rather than reconstructing the large incumbent? Compare residual-directed growth with reinsertion of removed incumbent atoms at matched size. |
| 3 | Reduce covariance rank C independently of B, with exact QR spectral compression and physical gauge alignment | At B719, C256 saves about 4,200 prep T: about 610M; C224 approaches 570M, before compiler/Q changes | Does the fixed integer dictionary have enough spectral rank margin? Exact full-signed PSD/rank-C lower bounds resolve coefficient-only continuation for this dictionary. |
| 4 | Dictionary-aware ANGLE: row permutation, shared parity features, nonlinear four-row lookup, and smaller J | Existing constructive J4 design gives T31608/Q454: 696.94M with banked lambda. A 1.25x win at Q454 requires total T <=25,693, hence ANGLE <=5,793 (about half the existing 11,708) | Do banked integer changes preserve compression relations, and can the bank shrink enough under the 128Q oracle cap? Compile exact table and run native oracle recovery/physical evaluation. |
| 5 | Spend unused one-body error on simpler one-body rows/relations, or jointly lower precision | One-body-only compression affects 76/795 rows: optimistic linear share of ANGLE is about 1,014 T, about 684.6M at fixed Q. Global brot7 also cuts Givens by 300 T and host thresholds may move | The host has independent 20 mHa bounds: the roughly 8 mHa unused margin cannot be transferred to two-body error. Test native one-body precision loss and resource sensitivity; simplification must produce an actual compiler saving. |
| 6 | Host-discontinuity continuation: coordinated B/C/R and J changes chosen by exact k2/address/prep thresholds, rather than monotone shape changes | B513->512->511 lowers Q478->462->460 as well as T; R2/B351/C175 can give about 488.4M conditional on physics; R4/B180/C88 about 395.7M | No host resource theorem proves the physical path exists. Map exact discrete jumps, then overlay measured continuation paths; this is a control strategy, not a new float shape sweep. |
| 7 | Preserve tensor but rotate physical factor gauge and trade error budget against Lambda | At existing resources universal normalization floor gives cost >682.185M, only 1.038x. Cannot reach 566.5M alone | Whether gauge optimization helps the smaller-table routes enough; test normalization after shape change, not as a standalone target claim. |

The top two are different directions through dictionary space: deletion from a
known feasible low-normalization state versus insertion of already-quantized
residual modes into a smaller dictionary. Neither repeats scratch float fits.
Rank and compiler probes supply independent controls and an attainable fallback.

Stopping: do not discard a route because its first draft loses. Resolve its
stated uncertainty on a declared finite proposal family, record exact-word
relaxation bounds, normalization tradeoffs and matched controls, and distinguish
that scoped negative evidence from global shape/architecture impossibility.
Any passing improvement receives full default native verification plus an
independent salt before local selection. No hash/salt search.

## Attribution and reproduction

Apache-2.0 native math/costs/oracle infrastructure: `from_host/verifiers/LICENSE.txt`.
Read-only reuse: chem-qfit exact QR projection and integer proposals; chem-audit
fresh signed target and physical derivatives; round-4 Codex signed Givens
encoder; chem-qrom structured ANGLE compilers. Source hashes will be recorded.
All scripts run from this directory after `source env.sh`.

## First numeric checkpoint

Native calibration reproduces banked error1 **11.995373 mHa**, error2
**19.989973 mHa**, Lambda **303.2035506**, gap **3.91504961**, lambda_eff
**48.56733747**. Fresh signed pair error agrees within 1.2e-16 Ha.

Exact unchanged-word PSD rank-C capacity at B719: C350 **19.516374 mHa**,
C340 **19.910161 mHa**, C320 **21.363025 mHa**, C256 **37.952099 mHa**,
C224 **55.803126 mHa**. Thus C320 or less is excluded for this dictionary;
this does not exclude different integer words. Conditional C224 resources
actually give **572.672M**, updating the rough Phase-1 estimate.

Sequential one-row deletion gives B703 rank-C error **20.182868 mHa**, B687
**21.347657 mHa**, B671 **23.943463 mHa**, B639 **56.586086 mHa**. The
unrestricted shared-span bound is already **20.284592 mHa** at B687.
Native physical-generator projection at B703 gives **20.627812 mHa** and
lambda_eff **48.98905**, whereas unconstrained error optimization has
lambda_eff **57.02470** after the damped gap projection. Accuracy alone is
not a cost win. Word repair is still required and is being tested.

Actual compiled J7 host resource screens give B575 **594.912M** (max lambda
for goal **46.24984**, below the universal bound), B543 **569.747M** (max
lambda **48.29266**), and B511 **523.516M** (max lambda **52.55732**).
These are conditional on banked normalization, not native passes.

Rounding one-body rows to an even-word grid gives error1 **123.108 mHa**;
rounding every row gives error2 **396.923 mHa**. The unused one-body margin
does not pay for a uniform lost precision bit in this basin. This only
rejects uniform rounding, not sensitivity-selected one-body simplification.

Commands: `"$CENSUS_PY" probe.py > results/probe.log 2>&1`, then
`bash run_searches.sh`. The second command runs bounded/unbounded word repair
from the same B703 physical-generator seed (60 iterations, seed 5, exact
rank-C solves), matched integer residual-growth/reinsertion at 479->511 and
639->671 (8 atoms per step), and fresh structured ANGLE compilation.

## Phase 2: first prototypes and follow-up

Bounded B703 repair accepts no word changes after line-limited exact weight
projection: **20.548637 mHa**, Lambda **319.9999983**. Three complete proposal
rounds reject 57 groups under Lambda/norm constraints. The unbounded matched
control accepts 60 rounds, improving **20.182868 ->20.057241 mHa**, but has
Lambda **398.65625**, lambda_eff **57.11693** before final offset optimization.
This first draft does not pass. It is being extended for 180 further rounds;
the initial endpoint is preserved, not overwritten.

At 639->671, residual growth reaches **34.984890 mHa**, while reinserting
incumbent atoms reaches **24.534255 mHa** under the same batch/rank-fit budget.
Their fixed-dictionary shared-span lower bounds are **34.352847** and
**23.618155 mHa**, respectively. These exclude coefficient-only rescue of
those dictionaries even without positivity/rank constraints. The growth
family does not rediscover the deleted incumbent geometry efficiently.

The initial 479->511 experiment used one-shot leverage pruning rather than
the sequential path, giving a poor seed (**1289.091 mHa**). Its endpoint is
preserved as seed-sensitive evidence. A follow-up uses sequential deletion
from the measured B511 state down to B479, then repeats both growth arms;
this avoids attributing one-shot seed damage to integer growth itself.

Fresh structured compilation of the banked words succeeds at **T31608/Q454**,
conditional cost **696,942,846.9115**, **1.016089x**, with unchanged physical
metrics. All nineteen oracles fit the 128Q cap. This is a compiler composition
of the banked physical fit with chem-qrom's architecture, not a new chemistry
basin or a 1.25x result. Full default native checks and an independent salt
are still required before selecting it.

Follow-up/reproduction: `bash run_followup.sh && bash run_physical.sh`.
All full native runs in the latter script use the shared verifier lock.

Exact host threshold detail: inner prep changes **17,318 ->16,967 T** on
B704->703, a 351T jump. The complete 16-row continuation additionally saves
ANGLE gates and small inverse-prep terms. The width cliff is split across
**B513/Q478 ->B512/Q462 ->B511/Q460**, not a single 18Q jump at B512.
`results/thresholds.json` retains every B479..719 point and each change.

The extended B703 repair is followed by a separate physical-generator
transfer and **6,000-step constrained weight/offset refit** using the exact
QR-compressed signed residual, error cap 19.990 mHa and strict gap target
3.91515. This tests whether a repaired dictionary's accuracy can actually be
retained at useful normalization, rather than judging raw unconstrained
weights alone. Any result still needs complete native evaluation.

Latest physical preflight: extended B703 reaches **19.982103 mHa** after
240 unrestricted repair rounds. The constrained feasible state has native
error2 **19.990000 mHa**, passes physical gates, but costs **807.506970M**
after finite offset optimization and exact host resource arithmetic. It is
not an improvement. Sequential-seed B511 growth gives **240.077776 mHa**
versus **244.992531 mHa** for matched incumbent reinsertion. The respective
span bounds are **240.021741** and **244.929077 mHa**, excluding those exact
word tables even with arbitrary real coefficients. Structured arbitrary-input
oracle checks pass every address in both directions under two independent
seeds. Full native scoring records will determine final selection.


<!-- audited-final-start -->
## Final measured result

**Verified improvement: 696,942,846.9115, 1.016089163x over the
banked 708,156,073.7876. The 1.25x target (566,524,859.0301) was not reached.**
Local selected candidate: `candidates/structured`; prior external selections
remain untouched. `selection.json` binds the exact candidate and both native
records. Native submission digest: `62f00bbfabeca29d4b9ae9f8f2fb69e9ed24f2e36db5883861438befd17a8216`.

| Candidate | Native status / shots | error1 mHa | error2 mHa | lambda_eff | T | Q | Product M |
|---|---|---:|---:|---:|---:|---:|---:|
| continuation703 | PASS / 49856 | 18.822246 | 19.990000 | 56.382923 | 29962 | 478 | 807.506970 |
| growth511 | REJECT / 37568 | 18.946003 | 240.077776 | 89.507770 | 23433 | 460 | 964.820361 |
| growth671 | REJECT / 47808 | 18.859718 | 34.984890 | 56.932096 | 28878 | 478 | 785.872660 |
| structured banked | PASS twice / 50880 each | 11.995373 | 19.989973 | 48.567337 | 31608 | 454 | 696.942847 |

Products for rejected rows are conditional arithmetic, not valid scores.
All rows use complete unchanged native ANGLE recovery, default 64 shots per
address, fixed physical shifts, native oracle semantics and every physical
gate. Full runs hold `/tmp/r6-fullverify.lock`. The selected candidate passes
two predeclared independent salts with identical digest and score. No salt
or hash search occurs.

## Top-route uncertainties and stopping decision

**Integer continuation:** the 240-round unrestricted repair reaches
19.982103 mHa at B703, demonstrating that changing integer
words can cross the first 16-row boundary. It does not establish the
low-normalization path to B511. The separate physical-generator transfer,
6,000-step constrained weight/offset fit, finite offset optimization and native
physical test quantify the normalization cost. A failed early 60-round
endpoint was extended rather than discarded. This family is not ruled out
globally; it remains an open route to deeper thresholds.

**Integer growth:** both residual-mode insertion and incumbent reinsertion
use exact integer atoms before all fitting, matched 8-row batches and rank-C
solves. A second 479->511 test repairs the poor one-shot seed by sequential
deletion first. Its endpoints are residual 240.077776 mHa (span bound 240.021741 mHa); reinsertion 244.992531 mHa (span bound 244.929077 mHa). At 639->671, residual growth
is 34.984890 mHa versus 24.534255 mHa for reinsertion. Independent span bounds
exclude all real coefficient-only rescue of these particular final growth
dictionaries. This resolves the finite candidate pool's capacity uncertainty
against it. Different seeds, refreshed pools, word moves and broader integer
growth remain unresolved, not globally excluded.

Independent CPU SciPy QR/eigensolves use exact native-decoded vectors and
the fresh signed target. Bounds have a 1e-8 Ha exclusion slack; these are
floating-point certificates, not interval-arithmetic proofs.

| Frozen research dictionary | Span bound mHa | PSD rank-C bound mHa | Saved residual mHa | Coefficient rescue excluded |
|---|---:|---:|---:|---|
| repair_B703_extended.npz | 18.836508 | 19.982103 | 19.982103 | no |
| growth511_sequential_residual_B511.npz | 240.021741 | 240.077776 | 240.077776 | yes |
| growth511_sequential_reinsertion_B511.npz | 244.929077 | 244.992531 | 244.992531 | yes |
| growth671_residual_B671.npz | 34.352847 | 34.984890 | 34.984890 | yes |
| growth671_reinsertion_B671.npz | 23.618155 | 24.534255 | 24.534255 | yes |

**Dictionary-independent rank audit:** every native two-body model is a PSD
matrix of rank at most R*C in the orthonormal symmetric pair coordinates.
The full fresh native spectrum gives necessary **R*C >= 290**
at 20 mHa, independent of dictionary, integer grid or factor gauge. This is
the PSD/rank approximation theorem evaluated in float64, not an interval
proof. It excludes standalone C256/C224 strategies even with different words
at R1. At B719 and this minimum C, the existing J7 compiler still gives
conditional **643.641068M** at banked normalization.
`results/global_rank.json` retains the adjacent threshold margins and full
spectral-tail numbers. This does not exclude combined B/rank/oracle changes.

**Other census routes:** the unchanged incumbent dictionary cannot reduce C
to 320 (21.363025 mHa optimum), far short of the C reduction needed for the
target. Uniform one-body precision reduction fails at 123.108293 mHa;
one-body and two-body bounds are independent, so unused one-body tolerance
cannot be reassigned to error2. The exact host map shows a split width cliff
at B513->512->511. The structured ANGLE fallback preserves the physical
factorization, increases total T but reduces persistent word storage by 24Q,
and improves the actual product by 1.583440%.
The compiler architecture is reused from chem-qrom and freshly fitted to the
banked words; this is not evidence for a new low-error physical basin.

These are finite, scoped candidate decisions. The census stops after the
declared continuation extension, constrained physical refit, matched growth
controls, fixed-dictionary capacity audits and complete native tests. It does
not claim that any shape, all integer growth, all ANGLE designs, or chemistry
as a whole is infeasible below the target. No experiment job is left running.

## Reproduction and artifact audit

From this directory:

```bash
source env.sh
"$CENSUS_PY" probe.py > results/probe.log 2>&1
bash run_searches.sh
bash run_followup.sh
bash run_physical.sh
"$CENSUS_PY" final_audit.py > results/audit.log
"$CENSUS_PY" report.py
```

`results/structured_components.json` checks all valid addresses, arbitrary
incoming words, both directions, clean scratch and phase under two independent
seeds. `results/audit.json` verifies all native digests/default workloads,
exact physical permutation, preserved banked NPZ and local file hashes.
`results/source_manifest.json` attributes and hashes read-only reused source;
`artifact_manifest.json` binds every numerical checkpoint and delivered
candidate. `NOTICE.md` attributes native Apache-2.0 and prior research code.
Fresh stochastic trajectories can differ numerically; delivered exact words,
weights, circuits and native digests bind the measured results.
<!-- audited-final-end -->
