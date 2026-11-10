# Round 5 integrated submission

The packaged submission passes fresh full-default native verification for all
three tasks. Its combined local score is **9.709873458706737**, up from the
round-four Codex archive's **9.367722153134604**, a **3.652449% improvement**.
Higher combined scores are better; lower individual task costs are better.

The archive is [submission.zip](submission.zip), with machine-readable costs
and provenance in [submission.summary.json](submission.summary.json).

| Component | Selected source | Verified shots | Q | Native task cost |
| --- | --- | ---: | ---: | ---: |
| Inverse square root | Round-four `zero_head_affine.kmb` | 200,540,928 | 90 | 182,956.5271606425 |
| DQI | Round-five `dqi-gap/gap_best_n70_m4095.kmb` | 102,400 | 1,004 | 1,435,930,717,449.9397 |
| Chemistry | Round-five `chem-audit/candidates/banked` | 50,880 | 478 | 708,156,073.7876049 |

The inverse circuit is unchanged from the round-four Codex selection. The DQI
decoder lowers its task cost by **10.153490%**, a **1.113009290743x** improvement.
The banked chemistry candidate supplies the independently verified small
normalization improvement from the second wave. Replacing only DQI while
retaining round-four chemistry would give a combined score of
**9.708085612889002**; the delivered archive includes both improvements.

## Native verification

All verifiers ran on files extracted from the delivered ZIP, using the unchanged
host code, default hash-derived seeds, default problem sizes and shot counts.
No seed override, shot reduction, or confidence adjustment was used.

- [Inverse verification](results/invsqrt_verify.json): PASS at all 200,540,928
  shots; worst error/tolerance 0.9796714395489131. Runtime: 160.4 seconds.
- [Chemistry verification](results/dfthcsos_verify.json): PASS at all 50,880
  shots. Errors are 11.995373 and 19.989973 mHa, with gap 3.915049606410321
  and effective normalization 48.56733747433665. Runtime: 13.7 seconds.
- [DQI verification](results/dqi_verify.json): PASS at all 102,400 shots, with
  **zero value errors and zero phase errors**. Mean executed Toffolis:
  1,195,802.1546875. Repetition factor: 1.1912604648115832, using the unchanged
  confidence parameter delta = 2^-128. Runtime: 142.9 seconds.

The fresh DQI score matches the candidate's original default-stream result.
Its prior default-plus-two-salt validation, adversarial checks, and residual
preservation report are copied under [results/source_evidence](results/source_evidence).
The decoder retains the documented rare degree-three and adjacent-double-gap
failures; its acceptance and score use the benchmark's approximate-decoder
confidence model. See the [source report](../dqi-gap/RESULTS.md) for that scope.

## Archive and score audit

The archive contains exactly 14 files: `invsqrt.kmb`, `dqi.kmb`,
`dfthcsos/submission.npz`, and `dfthcsos/angle_k_0.kmb` through
`dfthcsos/angle_k_10.kmb`. Its size is **160,459,941 bytes** (153.03 MiB).

Archive SHA256:

```text
d1b133764cd30750cbec60db0aecd1804c8779fc36907e9a2238cbc6d8d7f857
```

The [package manifest](results/package_manifest.json) records every source path,
member hash, byte count, source-verification report, and host verifier hash.
The [final archive audit](results/archive_audit.json) passed after reopening the
ZIP, reading and hashing every member, checking the extracted bytes, and
matching the native circuit hashes and canonical chemistry digest. Host
verifier and metric hashes remained unchanged.

DQI circuit SHA256:

```text
8a0675d51edd32a4eed4f5e35e2b3c5e13c83639bdb14142e36c119f7175fde9
```

The combined score uses the geometric mean implemented in
[the bundled metric](../../../from_host/metric.py), with the same local baseline
costs as the previous Codex archive:

```text
(2076774.233715 / inverse_cost
 * 24439463583672.93 / dqi_cost
 * 3355599588.315537 / chemistry_cost) ** (1/3)
```

The baseline improvement factors are 11.3511896293x for inverse square root,
17.0199462179x for DQI, and 4.7385028704x for chemistry. This is a local score;
Kaggle evaluates using its private split salts and associated baseline costs.

## Reproduction

From the repository root:

```bash
source initial_experiments/5/integration/env.sh
"$EXP_PYTHON" initial_experiments/5/integration/package_submission.py --workers 4
```

This repackages the selected bytes, extracts the ZIP, runs all three full
native verifiers sequentially under the shared full-verification lock, and
recomputes the summary and archive audit. Commands and timings are recorded in
[verification_runs.json](results/verification_runs.json).

To audit the existing archive and reports without rerunning simulation:

```bash
"$EXP_PYTHON" initial_experiments/5/integration/package_submission.py --audit-only
```

Original builders and reports remain in [dqi-gap](../dqi-gap/RESULTS.md),
[chem-audit](../chem-audit/RESULTS.md), and
[round-four Codex](../../4/codex/package_submission.py). Existing ignore rules
cover the archive, extracted `submission/` tree, logs, and generated results;
the packaging code and this report remain available to stage.
