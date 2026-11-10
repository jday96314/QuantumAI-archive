# Agent20 handoff report — 2026-10-25

**Current status.** The best packaged and audited submission scores **12.046515352810363**. Its independent-salt combination scores **12.04651477430945**. The target of 20 is unfinished. Milestones **10, 11 and 12** are saved; milestones 13–20 are missing. These are local scores under the supplied, unchanged host metric. No Kaggle upload was made.

A newer decoder finished both full native checks after the interruption. Combining its verified cost with the current other components gives **12.11607972269253** by default and **12.116078140634778** independently. **That combination has not been packaged or promoted.** The delivered `best/submission.zip` still contains the 12.046515 submission.

No experiment, GPU training job or verifier is running at this report's process check. The resumed goal manager reports `blocked`; no new technical or input dependency was identified. This report supersedes historical current-status claims in [CONTINUATION.md](CONTINUATION.md) and [RESULTS.md](RESULTS.md), including statements that the best is 10.128965.

**Saved submissions.** Milestones retain the first archive to cross each threshold.

| Artifact | Local default score | SHA-256 |
|---|---:|---|
| [best/submission.zip](best/submission.zip) | 12.046515352810363 | `f4aa6fc03c18c9028c7c2bb40954b0ffcd42defaa078e4a9789542cc02a6d78a` |
| [milestones/10/submission.zip](milestones/10/submission.zip) | 10.030435786719812 | `59bf461be714e0a4786fb0cc825d8b105622fec35ce193acd44ab1298004c2a0` |
| [milestones/11/submission.zip](milestones/11/submission.zip) | 11.747334147482889 | `c00ca03dfc2a5151204899e2d9f7b794640613697f7412f618174e4a056bf906` |
| [milestones/12/submission.zip](milestones/12/submission.zip) | 12.046515352810363 | `f4aa6fc03c18c9028c7c2bb40954b0ffcd42defaa078e4a9789542cc02a6d78a` |

The current best contains:

| Task | Selected source | Verified default cost | Full shots per run |
|---|---|---:|---:|
| Inverse square root | `public/phuc/invsqrt.kmb`, EXP225a | 138,095.07008539976 | 200,540,928 |
| DQI | `public/phuc/dqi.kmb`, EXP119a | 1,041,180,174,625.8771 | 102,400 |
| Chemistry | `candidates/paired_chemistry/` | 677,583,323.3862098 | 50,880 |

All three passed the unchanged full host verifier under default and predeclared independent salts. Both selected DQI runs had zero value and phase errors. The final archive passed the original metric readers and scoring audit. See [best/submission.summary.json](best/submission.summary.json), [best/results/archive_audit.json](best/results/archive_audit.json), [best/results/metric_preflight_and_score.json](best/results/metric_preflight_and_score.json), and [results/accepted_paired_chemistry_improvement.json](results/accepted_paired_chemistry_improvement.json).

**What worked: earlier local arithmetic.** Before the public-component improvements, the local combination reached **10.128965274112238**. An inverse-square-root row rewrite produced an 89-qubit circuit costing 173,447.66648. The repaired-H DQI endpoint used a Wronskian/Sylvester transform to reduce the point scan from 80 to 46 Toffolis. Its complete circuit used Q1004/static1123942 and cost 1.3557830990656e12 by default. Both full checks passed: one default phase error and zero independent phase errors.

Sources/evidence: `H_decoder.py`, `repaired_sylvester.py`, `fast_sylvester.py`, `results/repaired_H_build.json`, and `results/dqi_H_{default,independent}_verify.json`. The 10.128965 archive is preserved in `verified_before_public/`. This endpoint remains useful validated research but is no longer selected.

**What worked: finding and independently verifying better public components.** Public research found improvements that the local experiments had not achieved. I downloaded the published bytes, checked their hashes, retained attribution, and reran the unchanged host verifiers at full shots.

- [Ahmed Berat Özer's merged J-fraction decoder](https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression) shares one series reciprocal across two fraction levels. Its Q912 decoder passed twice locally, at default cost 1.071230662258e12. Combined with the then-published Phuc Thai EXP222a inverse circuit and our chemistry, it yielded **11.747334**, saving milestone 11. Sources and original LICENSE/NOTICE are in `public/ahmed_v3/`.
- [Phuc Thai's newer portfolio](https://www.kaggle.com/code/phucthaiv02/qcc-verified-portfolio-df-thc-retraining-study) supplied EXP225a inverse and EXP119a DQI bytes. EXP225a uses Q82/static1377/44 non-Clifford rotations. EXP119a uses Q898 with compact gap history. Both passed default and independent full checks locally. With the previous chemistry, they reached **11.933925**.

EXP225a SHA: `0f84f5bab3a02e6b357d4656553607174345a33d252cef643f9136db13212141`. EXP119a SHA: `c5da9cd4e2d73dd41417e539f2659009d646b3ee30ddb5ca494261378a502856`. Evidence: `results/public_phuc_*_verify.json`, `results/public_phuc_verification_runs.json`, and `parameters/phuc_release_pre_native.json`. That intermediate archive is in `verified_public_phuc/`.

**What worked: the paired chemistry lookup.** [paired_chemistry_angles.py](paired_chemistry_angles.py) adapts the paired bilinear-combine idea credited to Phuc Thai's fouraddr2 compiler to our existing compressed eight-bit bank:

```text
(a0 + d2)(a1 + d1) = a0*d1 + a1*d2 + d1*d2 + a0*a1   over GF(2).
```

Two lookup terms share one Toffoli. The known `d1*d2` correction is compiled into the base table, and `a0*a1` is shared. Physical NPZ bytes and loaded angle words are unchanged. The full block encoding saves **878 Toffolis**, **31,608 → 30,730**, with **454 qubits** and unchanged normalization. Chemistry cost falls **696,942,846.9115301 → 677,583,323.3862098**.

Every valid address passed arbitrary incoming output-word, repeated-application, scratch and phase checks under two seeds. Both full 50,880-shot native runs passed. This crossed 12. Evidence: `results/paired_chemistry_build.json`, `results/paired_chemistry_components.json`, `results/paired_chemistry_{default,independent}_verify.json`, and `results/paired_chemistry_verification_runs.json`. The delivered archive includes `public/PHUC_PAIRED_CHEM_NOTICE.md` as its NOTICE.

**Newest successful component, awaiting packaging: recovered gap flags.** [recovered_gap_decoder.py](recovered_gap_decoder.py) removes the persistent flag bank by recovering gap flags from fraction heads, and from the polynomial numerator during reverse conversion. It releases a source-derived temporary flag during scalar inversions. The full circuit is **Q886/static1100012**.

This is explicitly approximate. Its local rule uses `c_j == 0`, `a_j != 0`, and `c_(j-1) != 0`. Higher-order zero charts and certain consecutive ambiguous charts are unsupported. Guard HMR clears exceptional residual flag bits. It is **not an all-input literal-inverse construction**. Small paired checks passed on the tested physical and uniform syndromes, but deliberately forced unsupported zero charts failed. Those paired checks alone did not prove decoding truth; the subsequent full host checks did test the specified decoding phase.

| Full run | Value errors | Phase errors | Cost |
|---|---:|---:|---:|
| Default | 0 | 0 | 1,023,349,162,559.9081 |
| Independent | 0 | 0 | 1,023,349,270,684.8577 |

Each run used **102,400 shots**. Circuit: [candidates/recovered_gap_dqi.kmb](candidates/recovered_gap_dqi.kmb), SHA `a9222d2842aa6d53239ebbab35f449c4f4e864dd2b478151db73fa7b651646f1`. Evidence: `parameters/recovered_gap_pre_native.json`, `results/recovered_gap_dqi_build.json`, `results/recovered_gap_dqi_{default,independent}_verify.json`, `results/recovered_gap_verification_runs.json`, and `results/verify_recovered_gap.log`.

The build report's “requires native checks” status predates verification; the completed verifier JSONs supersede it. An early source adapter failed because of a missing `km` binding. That was corrected before the successful build. Historical failure text may remain in the raw build log; identify the successful result by its frozen source hashes, circuit hash and verifier reports. V1 source and native component bytes were archived separately.

**Decoder approaches that did not yield selected improvements.** Many reusable components passed native checks, but their whole-circuit costs or error rates lost. These are scoped results, not global impossibility claims.

| Direction | Outcome |
|---|---|
| Serial/blocked C- and J-fraction extraction | Small oracles and arithmetic identities passed; repeated long reciprocals and head reconstruction remained expensive. |
| Prefix/high-copy/mixed quotient decoders | A routed mixed design projected Q943/T967536, but a 100,000-vector census found 658 wrong decoded vectors. Complete costs and approximation penalties erased the gains. |
| Triple-gap head/frame/block | Canonical head passed at Q145/T2498 per direction. Integrated small blocks passed. Full extraction still cost **425579 T per direction**, before synthesis and scan. |
| Generalized degree-1/2 quotient reconstruction | Native tests passed, including arbitrary adjacent flags with skip priority. n35 v1 was Q989/T52454 per direction. Improved large v3 reconstruction was profiled around Q954/T48978, but this is reconstruction only. |
| Reversed-quotient cofactor / lookahead prefix | Components passed. A b6 N70 block cost Q1000/T102376, worse than the older generic prefix block's roughly 72k T. |
| Polynomial inverse / alternative phase endpoints | Useful primitives, but the complete extraction, conversion and scan budgets did not establish a better submission. |
| Field multiplier searches | Large finite searches did not improve the 42-Toffoli multiplier or the tested joint Forney quadratic decomposition. No optimality proof was obtained. |

See `triple_gap_*`, `generalized_gap_quotient.py`, `quotient_cofactor_frame.py`, `prefix_lookahead_quotient_block.py`, and associated `check_*`/`profile_*` files. Detailed inventories are in `results/head_quotient_research_manifest.json` and `results/generalized_quotient_research_manifest.json`. These are historical snapshots. Older recorder scripts assume an old best hash and must not be rerun unchanged against today's `best/`.

One interpretation trap: `serial_matches: 2` in the prefix-lookahead report compares different head encodings. It does not isolate tail decoding failure; that needs a supported-chart tail-only comparison.

**Chemistry approaches that did not pass.** Finite pruning/growth, physical refits and quantization-aware searches failed to improve the selected physical factorization at the required native accuracy.

- B703 pruning recovered accuracy, but worsened normalization and cost. Several deeper fixed dictionaries had projection floors excluding weight-only rescue.
- R4/B240/C88 positive fitting reached about **35.95 mHa continuous** after L-BFGS. Exact integer refinement reached **56.729335 mHa**, failing the **20 mHa** limit. The saved integer dictionary has an unrestricted union-span floor of **34.541746 mHa**: changing only weights cannot make those saved vectors pass.
- Fixed-factor ADMM, weighted ADMM, global atom replacement, and a weighted 20-inner-step variant did not improve that integer initializer. These runs are terminal; do not restart unchanged.
- R11/B64/C32 finished at **174.30 mHa continuous / about 267.22 mHa exact-word NNLS**. Its saved exact-word dictionary has an unrestricted span floor of **158.634857 mHa**. That excludes this dictionary, not the entire shape.
- The latest **40,000-step signed refinement** from the positive R4/B240 fit finished at **25.458428 mHa continuous**, **136.921887 mHa after exact-word refitting**. It failed and did not change the selected chemistry. Report: `results/signed_tensor_R4B240_from_positive.json`.

The earlier PSD generator factor is a useful matrix factor, not a compact angle dictionary. Necessary trace/rank constraints and fixed-dictionary bounds do not prove shape-wide infeasibility.

**Inverse-square-root approaches that lost.** Local fused-square, split-bank, direct signed-row and width/arithmetic variants were built and checked, but lost after lookup, arithmetic, phase routing, synthesis and live width were included. The rounded-square schedule, for example, used Q122/static2063/35 rotations and had an all-input numerical certificate, but cost more. The current winner is the independently verified public EXP225a artifact.

**Correction to old arithmetic notes.** Rechecked emitted counts: the retained full-word inverse costs **66 T**, Q45 including its 12 public source bits; the early-release version costs **102 T**, Q33 including the source. **34 T belongs to the specialized trace-inverse kernel**, not a full-word inverse. See `results/inverse_primitive_count_correction.json`. Native scores and actual emitter profiles remain valid; informal estimates using 34/66 for the two full-word variants need correction. An in-place field scale through `mon_slot` costs one 42-Toffoli multiplication plus measured erasure.

**Planned next steps.**

1. **Package the already verified recovered-gap decoder.** Check the frozen bindings in `parameters/recovered_gap_pre_native.json`. Prepare an accurate NOTICE for the modified decoder and its SomeJoint/Ahmed/manwithacat/anazemcev lineage; a Phuc-only decoder notice is incorrect for this replacement. Use `package_best.py` with the recovered decoder/default report, current Phuc inverse/default report, and paired chemistry/default report. Write to a separate directory, run `audit_best_metric.py`, and then `promote_verified.py` with the matching independent reports. This should deliver approximately **12.11608**. Full component verification is already complete; interruption alone does not justify repeating it.
2. **Try exact factoring of same-target inverse-square-root Toffoli groups.** This was the immediate unimplemented idea. `parameters/phuc_exp225.kmx` contains readable EXP225a gates. Its fused phase reads include `CX(s,t); CCX(a,s,t); CCX(b,s,t)`. The quadratic part factors to one `CCX(a+b,s,t)` using reversible CNOT basis changes and restoring the sources. General same-target quadratic forms may factor into fewer products of linear parities, without changing qubits, rotations or measurement behavior. **No optimizer for this idea was written or tested.** It needs an exact algebraic certificate and full native checks before selection.
3. Pursue larger architectural gains with complete cost accounting. From about 12.116 to 20 still requires approximately a **4.5-fold improvement in the product of task costs**. The merged extraction is a better starting point than repeatedly reconstructing full head frames. Chemistry needs new geometry or a substantially different feasible shape for a large gain; saved failed dictionaries cannot be rescued by more weight-only fitting.

**Environment and reproducibility.** The resumed session starts at the project root; work from agent20:

```bash
cd /media/shared/ExpansionDrive2/QuantumAI/initial_experiments/6/agent20
source env.sh
# EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
source chem_env.sh
# EXP_CHEM_PYTHON=/mnt/data01/miniconda3/envs/torch2.13/bin/python
"$EXP_CHEM_PYTHON" audit_best_metric.py best
```

Both interpreter paths existed at report time. Full verifiers use `/tmp/r6-fullverify.lock` and at most two workers. GPU experiments use `/tmp/r6-gpu.lock`. Keep new writes in agent20; prior experiments, host verifiers, reference data and metric stay read-only. Do not search salts or change sampling, tolerances, confidence formulas or costs. No subagents were spawned.

`package_best.py` accepts explicit task paths/results and NOTICE/LICENSE files, and saves crossed integer milestones. `promote_verified.py` preserves the previous best directory and checks the audit and independent report identities before promotion. `verify_frozen_release.py` runs unchanged host CLIs from a frozen manifest; archived driver versions preserve historical source hashes.

Milestone hashes were rechecked for this report. `results/report_snapshot.json` binds this report, the current best, milestones, newest decoder, and its source/verification evidence.
