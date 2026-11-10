# Experiment 4 / Codex

Final verified checkpoint: **9.3677221531**. Experiments were wrapped up at the
user's request; the **14.0** target has not been reached. No experiment jobs
remain running. [EXPERIMENT_REPORT.md](EXPERIMENT_REPORT.md) covers what worked,
what failed, verification evidence, and the remaining research leads.

[submission.zip](submission.zip) is the current archive.
[submission.summary.json](submission.summary.json) records the scores and digest;
[results/archive_audit.json](results/archive_audit.json) binds its fourteen
members to their verified bytes. The DQI member uses a new compact decoder.

| Task | Verified cost | Default shots | Evidence |
| --- | ---: | ---: | --- |
| Inverse square root | 182,956.527161 | 200,540,928 | [native result](results/zero_head_affine_verify.json) |
| DQI | 1,598,204,229,384.9265 | 102,400 | [native result](results/cubic_trace_product_n70_m4095_verify.json) |
| Chemistry | 708,547,388.955317 | 50,880 | [native result](results/norm_cycle2_signed_verify.json) |

The baseline archive scores 9.0096229084 and is retained under `baseline/`.
No host verifier or scoring code was modified.

The DQI circuit consumes syndromes into compact cofactors. Its new cubic trace
phase kernel uses88Toffolis per field point, down from92, and passes61440
exhaustive kernel shots. The complete decoder uses961qubits and1019865.8565mean
Toffolis, with zero value errors and380phase errors in102400native shots.
Its repetition factor is1.6968460733. The preceding9.3101204888archive is
preserved in`revisions/compact_product_slim/`. See
[CUBIC_TRACE_DQI_RESEARCH.md](CUBIC_TRACE_DQI_RESEARCH.md) for field identities,
source manifests, exhaustive checks and the rejected87Tvariant, and
[DEFERRED_COMPLETION_RESEARCH.md](DEFERRED_COMPLETION_RESEARCH.md) for the core.

The inverse circuit stores its constant coefficient after subtracting a cheap
affine function of the interval address. This reduces that coefficient from
27 to 23 bits. It restores the affine term before phasing, using a seven-bit
addition, a one-AND constant increment, and a parity-controlled CNOT. Two high
coefficient bits load late and clear early to provide arithmetic workspace.
The circuit retains 90 qubits and 34 synthesized rotations. Fixed-top-bit
restoration and measured cleanup of clean phase headroom reduce mean Toffolis
from1805.0004 to1782.0001. See
[ZERO_HEAD_AFFINE_RESEARCH.md](ZERO_HEAD_AFFINE_RESEARCH.md).
The preceding9.3326564887archive is preserved in`revisions/cubic_trace_product/`.
The previous 9.2156 archive is retained under `revisions/delayed/`, and the
9.1631 archive under `revisions/sign_anf/`.

The chemistry circuits use seven-word batches and two-entry lookups. For each
pair of rows they load the first row into the output and the XOR difference into
scratch, then apply the difference under the low address bit. Measurement and
an algebraic-normal-form phase repair erase the scratch. Six-bit lookup subtrees
keep each circuit within 128 qubits. Persistent angle storage falls from 104 to
56 qubits; total scored width falls from 526 to 478. An exact column-rotation
search further reduces normalization slightly. Its combined coefficient Gram
matrix is preserved to7.2e-13; native verification confirms the physical errors
and gap are preserved. See [RADICAL_SEARCH.md](RADICAL_SEARCH.md) for the ongoing
search beyond this small gain.

The selected family additionally uses an exact signed Givens gauge. Each
eight-bit angle is sign extended from seven independent bits, shrinking the
lookup difference bank from56 to49qubits. That permits larger unary subtrees
and lowers ANGLES cost to10604T, while preserving every physical metric.
The full native verifier and2048random-target checks per oracle pass.
See [ANGLE_GAUGE_RESEARCH.md](ANGLE_GAUGE_RESEARCH.md).

The selected dictionary further changes3505 integer angle words using coordinated
lattice proposals while retaining the factor-column gauge, then refits factors
and offsets. Its full native verifier, random-target oracle checks, exporter
rebuild and both ZIP audits pass. See
[NORM_DICTIONARY_RESEARCH.md](NORM_DICTIONARY_RESEARCH.md) for the exact selected
rebuild commands, frozen fit, unsuccessful alternatives and continuing search.
`package_submission.py` now selects this dictionary by default, and
`audit_submission.py` independently reopens and checks the resulting archive.

[check_selectswap.py](check_selectswap.py) also checks all eleven new XOR oracles
with random incoming target words, restored addresses, zero scratch, and correct
phase, beyond the host's streamed target states. Results are in
[results/selectswap_random_targets.json](results/selectswap_random_targets.json).

Rebuild and verify from this directory:

```bash
source env.sh
g++ -O3 compact_wires.cpp -o compact_wires
"$EXP_PYTHON" zero_head_affine_invsqrt.py
/mnt/data01/miniconda3/envs/torch2.13/bin/python ../../3/quantize_onebody.py \
  --fit parameters/norm_cycle2_final_fit.npz --out parameters/norm_cycle2_onebody.npz \
  --brot 8 --steps 3000 --target .012
"$EXP_PYTHON" ../../3/build_global_chem.py --fit parameters/norm_cycle2_onebody.npz \
  --out candidates/norm_cycle2_plain --batch 7
"$EXP_PYTHON" build_signed_gauge_chem.py --source candidates/norm_cycle2_plain/submission.npz \
  --out candidates/norm_cycle2_signed --lowbits 8
"$EXP_PYTHON" ../../../from_host/verifiers/invsqrt/verify_submission_invsqrt.py candidates/zero_head_affine.kmb --max-workers 4
"$EXP_PYTHON" ../../../from_host/verifiers/dfthcsos/verify.py candidates/norm_cycle2_signed
"$EXP_PYTHON" check_signed_gauge.py --candidate candidates/norm_cycle2_signed
"$EXP_PYTHON" build_cubic_compact_decoder.py
"$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py \
  candidates/cubic_trace_product_n70_m4095.kmb --max-workers 4
"$EXP_PYTHON" package_submission.py
"$EXP_PYTHON" audit_submission.py
```

The package command requires the saved full verifier results and checks their
hashes against the actual circuit files. Native simulation uses the bundled
kickmix extension and the existing Python 3.13 environment.

Failed experiments are retained as evidence, not included in the archive. The
smallest carry-accumulator schedule overflowed and failed at x=2. Single-parity
linear interpolation was explored but its gate-plus-rotation estimate did not
beat the selected quadratic circuit. The broad workspace search found no better
valid schedule than the inherited split sizes.

Further polynomial-arithmetic and dictionary-fit research is recorded in
[COMPACT_EUCLID_RESEARCH.md](COMPACT_EUCLID_RESEARCH.md),
[SERIES_RESEARCH.md](SERIES_RESEARCH.md), and
[LATTICE_RESEARCH.md](LATTICE_RESEARCH.md). Those experiments are excluded from
the archive.

The latest [blocked decoder experiments](BLOCKED_EUCLID_RESEARCH.md) include
quantum-tested compact quotient blocks and completion padding. Those earlier
profiles are superseded by the complete selected decoder documented in
[DEFERRED_COMPLETION_RESEARCH.md](DEFERRED_COMPLETION_RESEARCH.md). Ongoing
research is recorded in [RADICAL_SEARCH.md](RADICAL_SEARCH.md).

[NORMALIZATION_BOUND.md](NORMALIZATION_BOUND.md) derives the global native
normalization lowerbound46.7862, directing larger gains toward gates and width.

[GROUPED_FACTOR_RESEARCH.md](GROUPED_FACTOR_RESEARCH.md) records alternative
chemistry shapes whose quantized dictionaries did not reach the physical error
limit; the R28floating-point fit reaches17.02mHa before rounding. A new direct
Boolean phase model in`boolean_polynomial_phase.py` matches the selected integer
inverse-square-root model on32857sampled/boundary inputs. Its795AND-node graph
has not been scheduled into a quantum circuit; it establishes no resource gain.

The inherited decoder adapts Apache-2.0 code from
[andigu/uncompute](https://github.com/andigu/uncompute), revision
`64b5a13aa249f1360846657392e6c66817b26d0f`. Its source and license remain in
`../../3/vendor/uncompute`. See [CONTINUATION.md](CONTINUATION.md) for remaining
research and exact checkpoint details.

[PHASE_AND_GROUPED_RESEARCH.md](PHASE_AND_GROUPED_RESEARCH.md) records the
rejected cubic arithmetic circuit, Goldschmidt models, Walsh screen and grouped
quantization probes. New decoder look-ahead and series-block experiments are
recorded in[LOOKAHEAD_RESEARCH.md](LOOKAHEAD_RESEARCH.md).
