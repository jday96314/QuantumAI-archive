ACCEPTED: trimmed decoder

Accepted on lower T at equal Q and exact paired-seed decoding, not on failure-count luck. Full native Q=961, static T=994358 (25596 fewer than incumbent). The schedule and cubic scan are unchanged.

Builder: `build_accepted_decoder.py`; KMB: `candidates/trimmed_decoder.kmb`; SHA-256 `91f3e6680ac51af5b56b8baf8b051b2411091d5a3fc845fd53bface0ee44582d`.

| Full 102400-shot run | Value / phase failures | R | Mean T | Cost |
| --- | ---: | ---: | ---: | ---: |
| Incumbent, paired seed 2026102405 | 0 / 395 | 1.71157557352 | 1019866.334296875051 | 1.61207822277e+12 |
| Trimmed, paired seed 2026102405 | 0 / 395 | 1.71157557352 | 994270.334296875051 | 1.57161924025e+12 |
| Trimmed, default | 0 / 382 | 1.69881437737 | 994269.954453125014 | 1.55990093713e+12 |
| Trimmed, independent salt | 0 / 379 | 1.69586139626 | 994270.035234375042 | 1.55718955022e+12 |

All four runs used the unchanged full-default host verifier under `flock /tmp/r6-fullverify.lock`, with at most two CPU workers. The paired runs have identical failure counts, R and Q, and exactly 25596 fewer mean T. JSONs: `results/incumbent_seed_2026102405.json`, `results/trimmed_seed_2026102405.json`, `results/trimmed_default_verify.json`, `results/trimmed_independent_verify.json`; acceptance: `results/accepted_decoder.json`. Independent salt: `dqi-core-trimmed-independent-20261024-a`.

Rebuild from repository root: `source initial_experiments/5/dqi-core/env.sh`, then `"$EXP_PYTHON" initial_experiments/5/dqi-core/build_accepted_decoder.py`. Re-run verification after the incumbent paired report exists with `"$EXP_PYTHON" initial_experiments/5/dqi-core/verify_and_accept.py`. Commands and fixed sampling plan are recorded in `parameters/verification_plan.json` and `results/*_command.json`.

The earlier research-only stopping result below is superseded by this full-circuit acceptance. The protected round-4 incumbent remains unchanged. Work stops here.

---

# Compact Euclid core research

The round-4 selected artifacts are read-only and remain selected. Reference
overall score 9.3677221531; inverse 182956.5272, DQI 1598204229385,
chemistry 708547388.9553. No verifier, confidence calculation or scoring changed.

Starting evidence: selected core 649802 T, scan 370152 T; full selected width
961 Q, not 852 Q. The compact input/state itself is 840 Q. Target: reduce core
at unchanged selected peak, and report separately against the tighter 852 Q
target. Counts and estimates will be distinguished; no replacement is accepted
without full native 102400-shot verification and independent salts.

Hypotheses tested: dense triangular tail transforms dominate; selective
measurement erasure can only help if a forward temporary remains available
at the reverse boundary; larger early blocks may amortize tail work but can
exceed the head workspace limit.

## Outcome and first numeric checkpoint

**Stopped under the <10%-reducible scope clause.** Best counted core:
**624206 T = 312103 forward + 312103 reverse**, saving **25596 T (3.939%)**,
at unchanged selected **961-Q** peak. The ≥20% target and 852-Q target were
not achieved. No full decoder was emitted, accepted or promoted. This result
does not bound other decoder designs. The complete analysis and algorithmic
options are in [DESIGN_MEMO.md](DESIGN_MEMO.md).

The first checkpoint reproduced 324901 T per direction. Series/middle products
cost 162582 T per direction (50.04%). Final exclusive inverse attribution is
16320 T; the earlier scratch count included work inside inverse contexts,
which is now separated correctly. No total gate count changed in that correction.

| Experiment | Forward T | Reverse T | Peak Q | Status |
| --- | ---: | ---: | ---: | --- |
| Selected exact replay | 324901 | 324901 | 961 | Reference reproduced |
| Trimmed products + budgeted inverse retention | 312103 | 312103 | 961 | Native components pass |
| Optimal tested width-feasible schedule | 312103 | 312103 | 961 | Same block boundaries |
| Original block4, serial-after-21 | 307837 | Not separately counted | 984 | Width failure |
| Optimized standalone early block4 | 20714 | Not counted | 968 | Width failure; partial component |

Changes: trim zero rows/columns before choosing middle-product algorithm,
including offset suffix maps (11214 T saved per direction); retain inverse
chain in 44 eligible contexts per direction (1584 T saved per direction).
No forward temporary history survives at the scan boundary for an additional
mechanical asymmetric reverse erasure; estimated additional saving is 0 T.

Checks: **28 native circuits × 512 shots = 14336 component shots**, independent
polynomial values, inverse/composition, q=1/2/4 middle maps, completion cases,
source/flag preservation, clean scratch and zero phase. Twelve native block
counts independently equal the counter's results at k=1,19,22,34.
The count screen covers **120 local alternatives, both directions**, followed
by dynamic programming under Q≤961. It returns block3 at k=1,4,...,19 and
serial at k=22..34. Native artifacts and SHA-256 digests are retained locally.

Failures: a full running-ratio recurrence costs 175953 / 145758 / 81660 T
for three steps starting at k=1 / 10 / 19, versus 15687 / 32769 / 35007 T
for the optimized head-block design. Its lower width does not offset repeated
whole-series inversions. This rejects that recurrence only. No in-place
(sigma,omega) adapter was produced; its break-even budget for an 80-T scan
is 16380 T per direction, including completion handling.

**Measured versus estimated:** attribution/schedule values are exact native
emitter counts; component KMB counts and simulated values are native checked.
Full decoder counts are projected by composition. Static total falls from
1019954 to **994358 T** including the unchanged 370152-T scan (2.5095%).
Assuming R=1.2 and approximating mean T by static T gives **1.1019725934e12**
for R Q² T at Q=961. R=1.2 is assumed; selected measured R is 1.6968460733.
At hypothetical Q=852 the projection is 8.6617013956e11, but that width is
not achieved. Full 102400-shot verification and independent salts were not
run because this research stops below the requested reduction threshold.

Machine-readable checkpoint: [results/checkpoint.json](results/checkpoint.json).
Detailed attribution, schedule and native reports are in `results/`;
the reviewable native prototypes are in `candidates/`.

Reproduction (from the repository root):

```bash
source initial_experiments/5/dqi-core/env.sh
"$EXP_PYTHON" initial_experiments/5/dqi-core/attribute_core.py
"$EXP_PYTHON" initial_experiments/5/dqi-core/attribute_core.py --trimmed --retain
"$EXP_PYTHON" initial_experiments/5/dqi-core/attribute_core.py --block 4 --serial-after 21 --forward-only
"$EXP_PYTHON" initial_experiments/5/dqi-core/native_components.py middle --shots 512
"$EXP_PYTHON" initial_experiments/5/dqi-core/native_components.py blocks --shots 512
"$EXP_PYTHON" initial_experiments/5/dqi-core/schedule_screen.py
"$EXP_PYTHON" initial_experiments/5/dqi-core/ratio_probe.py
"$EXP_PYTHON" initial_experiments/5/dqi-core/summarize.py
```

Run sequentially or at most two processes concurrently; environment caps BLAS
to one thread. No full verifier is invoked. Any future full-default verification
must use `flock /tmp/r6-fullverify.lock` and native default shots plus independent
salts before acceptance. All output paths stay inside this directory.

Attribution reuses the round-5 `dqi-hybrid/profile_selected.py` replay structure,
`dqi-phase/baseline_audit.py` selected-candidate identification, and round-4
`profile_compact_parts.py` product-term counting. Arithmetic is reused via
`import_published.load()` from andigu/uncompute, Apache-2.0, revision
64b5a13aa249f1360846657392e6c66817b26d0f. Original source/license remain
in the read-only round-3 vendor tree. No upstream source is copied or modified.
