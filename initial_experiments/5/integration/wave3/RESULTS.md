# Round 5 wave-three integrated submission

**Verified local score: 9.811059643318**, up from **9.709873458707** (1.042096% higher).
The deliverable is [submission.zip](submission.zip). The previously submitted archive and its report remain unchanged one directory above.

## Integrated improvements

| Task | Previous task cost | New native task cost | Cost reduction |
| --- | ---: | ---: | ---: |
| invsqrt | 182,956.5271606425 | 182,956.5271606425 | 0.000000% |
| dqi | 1,435,930,717,449.9396972656 | 1,414,354,472,410.0034179688 | 1.502597% |
| dfthcsos | 708,156,073.7876049280 | 696,942,846.9115301371 | 1.583440% |

The accepted wave-three middle-product emitter removes **23,352 static Toffolis** from the existing gap-aware DQI decoder, at the same **1,004-qubit peak**. The gap-aware schedule, dispatch, original low-workspace inverse, and phase scan are retained. The complete compiled count is **1,172,538 static Toffolis** (801,388 core + 371,150 scan).

Chemistry uses `chem-census/candidates/structured`: the content-aware QROM compiler composed with the banked physical fit. It uses 19 ANGLE circuits, 31,608 host Toffolis and 454 host qubits. This already incorporates the overlapping QROM and normalization gains.

The inverse remains the round-four `zero_head_affine.kmb`; neither inverse report supplies a cheaper accepted replacement.

## Review of all ten reports

| Report | Integration decision |
| --- | --- |
| [chem_audit_pruning_and_submission_report](../../brainstorming/results/wave3/claude/chem_audit_pruning_and_submission_report.md) | Zero subtree-pruning savings; its banked physical improvement is already present. Its packaged trimmed DQI is slower than the accepted gap decoder. |
| [chem_census_factorization_and_compilation_report](../../brainstorming/results/wave3/claude/chem_census_factorization_and_compilation_report.md) | Selected: verified structured compiler applied to the banked physical fit, 696942846.9115301 chemistry cost. |
| [dqi_core_trimmed_decoder_report](../../brainstorming/results/wave3/claude/dqi_core_trimmed_decoder_report.md) | Ported only the accepted trimmed middle-product emitter to the gap decoder; directly replacing the gap decoder with the old trimmed architecture would regress cost. |
| [dqi_dispatch_exception_correction_report](../../brainstorming/results/wave3/claude/dqi_dispatch_exception_correction_report.md) | No complete accepted replacement; resource projections omit required work. |
| [dqi_census_decoder_architecture_report](../../brainstorming/results/wave3/claude/dqi_census_decoder_architecture_report.md) | Working research components and conditional budgets; no full verified replacement. |
| [inv_census_inverse_architecture_report](../../brainstorming/results/wave3/claude/inv_census_inverse_architecture_report.md) | No lower-cost verified inverse; retain the current circuit. |
| [chem_qrom_content_aware_compilation_report](../../brainstorming/results/wave3/gpt/chem_qrom_content_aware_compilation_report.md) | Verified 697327965.8381622 compiler gain; already composed with better banked normalization in chem-census, so do not count it twice. |
| [chem_tree_topology_and_quantization_report](../../brainstorming/results/wave3/gpt/chem_tree_topology_and_quantization_report.md) | No physically valid cheaper replacement. |
| [dqi_core_phase_aware_decoder_report](../../brainstorming/results/wave3/gpt/dqi_core_phase_aware_decoder_report.md) | Incomplete gap handling and conversion in candidate architectures; no accepted full-task gain. |
| [inv_phase_inverse_square_root_report](../../brainstorming/results/wave3/gpt/inv_phase_inverse_square_root_report.md) | Both full native candidates pass but cost more than the current inverse. |

## Verification and score

The DQI candidate was fixed from resource counts before any native sampling. Its source manifest records the accepted gap sources and the exact emitter from the accepted trimmed-decoder build. No budgeted inverse-retention change was ported.

| DQI run | Shots | Value / phase errors | Mean executed Toffolis | Native cost |
| --- | ---: | ---: | ---: | ---: |
| default | 102,400 | 0 / 1 | 1,172,449.431093750 | 1,414,354,472,410.003418 |
| independent | 102,400 | 0 / 1 | 1,172,450.077343750 | 1,414,355,251,997.210938 |
| paired | 102,400 | 0 / 0 | 1,172,450.154687500 | 1,407,889,411,467.635498 |

The paired run uses the original gap decoder's default hash-derived seed. It reproduces the same decoding counters, repetition factor, and width, with exactly **23,352 fewer mean Toffolis**. Default and independent-salt runs pass separately; the default supplies the score. Confidence remains delta = 2^-128 and each full run uses 102,400 shots.

The new default stream has **1 phase error and 0 value errors**. A further full 102,400-shot replay of the original gap decoder on that same stream reproduces those counters and the repetition factor, with exactly 23,352 more mean Toffolis. The reported score uses the default stream with its actual confidence penalty. See the [controlled stream comparison](results/new_default_stream_comparison.json).

The full-mask adversarial test passes 1,024 cases covering all gap stages, immediate completion, zero/one-error, weight-35, and random controls, including twice-executed restoration. All five known residual examples pass 640 independent-mask input/scratch and twice-execution restoration checks, retaining the original phase-error counts. The decoder still approximates the rare degree-three and adjacent-double-gap cases.

Fresh chemistry verification passes all **50,880 shots** at default settings, matching the two independent native passes in its source experiment. The inverse's unchanged bytes reuse their **200,540,928-shot PASS** evidence. No reduced workloads or favorable-salt selection are used.

The package audit reopens the ZIP, reads every member to EOF, verifies hashes against the selected source files and extracted artifacts, and matches the native DQI/inverse circuit hashes and canonical chemistry digest. The unmodified bundled metric's circuit and NPZ readers pass, and its `score` function reproduces the combined score from the hash-bound native costs. This is a local default-stream score using the same baseline costs as the prior integration.

Archive size: **154,589,243 bytes**. Archive SHA256:

```text
b521c86dfdf830ada6e0f9cbfc88a5e58c762be3a2feabaa50c0479324652fa1
```

Machine-readable records:

- [Submission summary](submission.summary.json)
- [Native DQI acceptance](results/gap_trimmed_acceptance.json)
- [Adversarial restoration](results/gap_trimmed_adversarial.json) and [known residuals](results/gap_trimmed_residual.json)
- [Archive audit](results/archive_audit.json) and [unmodified metric readers/score](results/metric_preflight_and_score.json)
- [All report decisions and source hashes](results/report_selection_review.json)

## Reproduction

From the repository root:

```bash
source initial_experiments/5/integration/env.sh
"$EXP_PYTHON" initial_experiments/5/integration/probe_gap_trimming.py
"$EXP_PYTHON" initial_experiments/5/integration/build_wave3_decoder.py
"$EXP_PYTHON" initial_experiments/5/integration/validate_wave3_decoder.py
"$EXP_PYTHON" initial_experiments/5/integration/check_wave3_default_stream.py
flock /tmp/r6-fullverify.lock "$EXP_PYTHON" from_host/verifiers/dfthcsos/verify.py initial_experiments/5/chem-census/candidates/structured > initial_experiments/5/integration/wave3/results/chemistry_verify_default.json
bash initial_experiments/5/integration/package_wave3.sh
```

`package_wave3.sh` packages the accepted native-verified inputs and audits the exact archive. It reuses the already completed full native checks after binding their digests to the packaged bytes. Re-running that script alone does not repeat simulation. Reused field arithmetic retains its upstream Apache-2.0 attribution in the original experiment reports.
