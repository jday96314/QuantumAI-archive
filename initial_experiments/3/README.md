# Experiment 3

The goal is reached: the fully verified overall score is **9.009622908**.
All three unchanged host verifiers report `PASS` at their default problem sizes
and shot counts.

- Ready-to-upload archive: [submission.zip](submission.zip).
- Exact scores, source results, and SHA-256: [submission.summary.json](submission.summary.json).
- Milestones 5 through 9 are preserved in separate directories under `milestones/`.
- Milestone 9 retains the original checkpoint; the current archive matches `revisions/mixing_mu1e7/submission.zip`.
- [Experiment report](EXPERIMENT_REPORT.md): successful and unsuccessful ideas for all tasks, plus the matched sign-mixing study.

| Task | Verified cost (lower is better) | Improvement over baseline | Verified shots |
| --- | ---: | ---: | ---: |
| Inverse square root | 190,711.955455 | 10.889586× | 200,540,928 |
| DQI | 1,612,737,419,463.0198 | 15.154025× | 102,400 |
| DFTHC-SOS chemistry | 757,163,489.498628 | 4.431803× | 50,880 |

Archive SHA-256: `6ebc780a1c8027f26080e47818ea5db63dde58ce9ffbf90f46c4e7411b904800`.

## Follow-up mixing experiment

A weak error-scaled mixing penalty (μ=10⁻⁷) slightly improved the matched
unpenalized refit. The original μ≈10⁻⁶ suggestion did not produce the predicted
normalization reduction. The current chemistry cost is 0.036031% below the old
submission, and 0.001154% below the matched control. Coefficient conventions,
full native results, and limitations are in the [report](EXPERIMENT_REPORT.md)
and [experiment summary](results/mixing_experiment_summary.json).

## Final changes

The inverse-square-root circuit uses 94 qubits. Hybrid adders trade a few gates
for lower workspace use. An unused curvature coefficient is measured during
Horner evaluation; its measurement phase is repaired when that coefficient is
reconstructed for the inverse evaluation. Idle output bits supply integration
carries, and private-wire lifetime compaction removes allocator gaps. The
coefficient table and the computed phase function remain fixed.

The chemistry factorization has `(R, B, C) = (1, 719, 350)` and eight-bit angle
words. Refining its factors under the physical error and SOS-gap constraints
lowers the verified effective normalization to 48.8951873178. Six angle-loading
circuits use batches of 13 words. The verifier reports 29,440 Toffolis and 526
scored qubits. One-body angle optimization is included in the verified files.

The DQI circuit consumes syndrome storage as it becomes available and contracts
inverse-parity operations through subfields. Its verified repetition factor is
included in its score. Its full result is linked below.

## Verification and reproduction

The final archive audit checks that all nine archive members match their
verified source files, the chemistry submission digest matches the verifier,
the archive checksum matches the summary, and milestones 5–9 retain their
recorded checksums: [final_archive_audit.json](results/final_archive_audit.json).

The complete host results are:

- [Inverse square root](results/compactdefer_k4_q94_verify.json).
- [DQI](results/stream10_offset3_select_noscale_compact_verify.json).
- [Chemistry](results/mixing_mu1e7_verify.json).

Run from the repository root after sourcing `initial_experiments/3/env.sh`.
`build_compact_invsqrt.py --output PATH.kmb` reproduces the final inverse circuit;
its bytes were checked against the verified candidate. `build_global_chem.py`
rebuilds the chemistry circuit family from
`candidates/mixing_mu1e7_onebody.npz`. `package.py` computes the geometric
mean of the verified baseline improvements and rejects a milestone archive
unless its threshold is reached. The research history and unsuccessful branches
are retained in [CONTINUATION.md](CONTINUATION.md).

## Published decoder adaptation

`vendor/uncompute` is the Apache-2.0 upstream repository
<https://github.com/andigu/uncompute>, revision
`64b5a13aa249f1360846657392e6c66817b26d0f`; its license is retained there.
`import_published.py`, `published_dqi.py`, `translate_published.cpp`, and
`build_published.py` adapt its Reed–Solomon decoder, replacing the fixed-mask
sweep with a runtime-mask circuit and converting between field bases with
reversible linear maps. The final adaptation passes the unchanged host verifier.
The upstream polynomial reciprocal derivative construction is documented at
<https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html>.
