# Experiment report

Experiments were wrapped up at the user's request. The best audited submission
scores **9.367722153134604**. The requested **14.0** score was **not achieved**.
No experiment jobs remain running.

The campaign improved the inherited local checkpoint from 9.00962290835304 to
9.367722153134604, a **3.97% score increase**. Several substantial circuit
changes worked, but their practical gains were much smaller than their raw
gate or qubit reductions suggested. The most useful outcome beyond the
submission is a tested collection of reversible arithmetic components, plus
evidence identifying which proposed architectures fail accuracy, cost, or
reversibility requirements.

**The final deliverable.** [submission.zip](submission.zip) contains the selected
inverse, DQI decoder, and chemistry submission. Its SHA256 is
`9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68`.
[submission.summary.json](submission.summary.json) records the task costs and
score. Both the [packaging audit](results/archive_audit.json) and
[independent archive audit](results/archive_independent_audit.json) pass.
They bind all 14 members to the native verification results and canonical
chemistry parameter digest. The previous 9.3326564887 archive and its audits are
preserved in [revisions/cubic_trace_product/](revisions/cubic_trace_product/).
Earlier checkpoints remain under `revisions/`; the inherited checkpoint is in
`baseline/`.

| Task | Inherited local cost | Final cost | Cost reduction |
| --- | ---: | ---: | ---: |
| Inverse square root |190,711.955455|182,956.527161|4.07%|
| DQI decoder |1.612737419463e12|1.598204229385e12|0.90%|
| Chemistry |757,163,489.498628|708,547,388.955317|6.42%|

Costs are compared within each task. The overall competition score is the
geometric mean of the three improvements over the competition references.

**What worked: inverse square root.** The selected circuit is
[zero_head_affine.kmb](candidates/zero_head_affine.kmb). It passes all
**200,540,928 native shots**, uses 90 qubits and 34 rotations, and has
1782.000090874218 mean Toffolis, 1793 static Toffolis. Its maximum
error/tolerance ratio is 0.9796714395489131 at x=958102.

The main storage improvement subtracts a cheap affine function of the interval
address from the lookup's constant coefficient. That coefficient then fits in
23 bits instead of 27. Delaying high coefficient bits and releasing curvature
bits early provides arithmetic workspace without widening the circuit.

Two final changes reduced the preceding circuit's 1816 static Toffolis to 1793.
An exhaustive integer-model check established that the normalized phase
mantissa always has bit 26 equal to one. Restoring the affine offset modulo
2^26 therefore preserves the correct result while shortening the addition.
The constant restoration uses a one-AND increment and the parity contribution
uses a CNOT. Output shifts then exploit their clean private headroom: the
forward move computes an AND, while its reverse uses measurement and a
Clifford phase repair. This avoids paying for the reverse controlled swap.

The complete selected inverse was rebuilt independently and produced identical
KMB bytes. See [native verification](results/zero_head_affine_verify.json),
[rebuild evidence](results/zero_head_affine_rebuild_check.json),
[source manifest](parameters/zero_head_affine_source_manifest.json), and
[the derivation and rebuild instructions](ZERO_HEAD_AFFINE_RESEARCH.md).

**What worked: DQI.** The selected circuit is
[cubic_trace_product_n70_m4095.kmb](candidates/cubic_trace_product_n70_m4095.kmb).
It consumes the 70 syndrome words into compact polynomial cofactors, reducing
the total width to 961 qubits. Its core uses blocked Euclidean updates, switches
to serial updates late in the computation, and postpones completion handling
so that it can reuse syndrome storage.

The accepted phase kernel represents GF(4096) as a cubic extension of GF(16).
A trace identity reduces inverse evaluation to a small norm and a five-AND
GF(16) inverse. Runtime classical masks are incorporated through reversible
linear maps. This reduces the scan kernel from 92 to 88 Toffolis per field
position. The kernel passes 61,440 checks, including all 4,096 divisors under
12 basis masks, additional masks, zero cases, input restoration, scratch, and
phase. The complete decoder was
also independently rebuilt with identical bytes.

The full native decoder passes 102,400 shots with **zero value errors and 380
phase errors**. Those phase errors are permitted and charged by the native
scoring rule; this is not an exact decoder claim. The resulting repetition
factor is 1.6968460732519504 and mean Toffolis are 1,019,865.856484375.
The inherited decoder used more qubits but fewer gates and a lower repetition
factor. Most of the width benefit was consumed by those increases, leaving
only a 0.90% task-cost improvement.

See [native verification](results/cubic_trace_product_n70_m4095_verify.json),
[rebuild evidence](results/cubic_trace_product_rebuild_check.json),
[field-kernel research](CUBIC_TRACE_DQI_RESEARCH.md), and
[compact-core research](DEFERRED_COMPLETION_RESEARCH.md).

**What worked: chemistry.** The selected directory is
[norm_cycle2_signed](candidates/norm_cycle2_signed/). Its 11 angle oracles use
seven-word batches. A signed angle gauge preserves the physical projectors
while making two output bits dependent, reducing the independently loaded
payload. That permits larger lookup subtrees within the 128-qubit oracle cap.
Each selected oracle uses 123 qubits; the complete scored preparation uses 478.

Normalization-aware integer dictionary updates and factor/offset refits also
improved the physical model. The final submission passes 50,880 native shots
with 30,504 Toffolis and effective normalization 48.59417497035283. One-body
error is 0.011986315442796802 Ha; two-body error is 0.0199953044905532 Ha.
The two-body error is close to the 0.020 Ha limit, so its frozen quantized
parameters matter. The SOS gap is 3.915157472365536 Ha, also close to the native
floor of 3.914898 Ha.

All 11 oracles additionally pass 2,048 random incoming-target shots each, including
address restoration, zero scratch and correct phase. An independent exporter
rebuild reproduced the oracle bytes and parameter arrays. See
[native verification](results/norm_cycle2_signed_verify.json),
[signed-gauge research](ANGLE_GAUGE_RESEARCH.md), and
[normalization/dictionary research](NORM_DICTIONARY_RESEARCH.md).

**What did not work.** The following outcomes are specific to the tested
implementations and searches. A failed fit or bounded search is not a general
impossibility proof.

| Approach | Observed outcome | Why it was excluded |
| --- | --- | --- |
| Direct cubic inverse approximation | Corrected circuit passed 131,072 independent quantum checks; 94Q/2458T. | More expensive than the selected quadratic circuit. Earlier versions also exposed accumulator overflow, which was corrected before the passing result. |
| Goldschmidt seed/residual correction | Whole-domain numerical models now pass using actual Booth-row rounding. A reversible split-square component passes 16,384 cases in each direction and both composition orders. | No complete competitive circuit was established. The smaller square component alone costs 425T per direction, before products, lookup replay, normalization and phase. |
| Direct Boolean phase compilation and measured pebbling | The 795-AND arithmetic model matches the integer reference. Generic quantum pebbling helpers pass their tests. | The tested schedules required excessive workspace or recomputation. ABC resynthesis increased the AND count. |
| Walsh phase expansion | 2,048 terms meet the numerical tolerance; 1,024 do not. | The required rotation count is too large before adding exponent and parity handling. |
| An 87T DQI scan kernel | Full 102,400-shot native verification passed, with lower gate count. | Its scored cost was 1.61569700396598e12, worse than the selected 88T variant. Circuit-dependent sampling changed the measured phase-error penalty. This is not evidence that the cheaper kernel is intrinsically less accurate. |
| Full-series and terminal Euclid blocks | Reversible components pass. All 78 tested late-block schedules lose to the current serial schedule. | Lower component workspace did not offset greater gate count or reduce the complete decoder's earlier width peak. |
| Six-moment look-ahead coordinates | 16,384 cases per direction pass, including forced zero pivots, singular extensions, scratch and phase. An opening-block integration also passes. | The tested initial four-degree prefix costs 28,114T/Q970 versus 20,360T/Q943. A general global-cofactor adapter remains unbuilt. |
| Thin cofactor updates | Smaller preliminary resource counts. | An arbitrary-input round trip failed at shot 199. These counts are invalid as submission evidence. |
| Smaller chemistry dictionaries | Several floating-point fits improve substantially. | Low-precision quantization and refitting missed the 20 mHa physical error gate. Higher precision increased the native lookup workspace enough to lose the cost advantage. |
| Grouped chemistry factors | R28/B64/C32 reaches 17.0208 mHa before quantization. | Direct 7-bit quantization gives 391.182 mHa; even a larger fixed-dictionary covariance refit remains 229.739 mHa. The tested 24-group fit misses the error gate even before quantization. |
| Extra atoms and sparse coordinate dictionaries | Some fits or normalization values improve. | Added lookup resources outweigh the gain, or the quantized residual remains invalid. The full one-hot/pair dictionary also has a poor normalization/resource certificate. |
| Cheaper field multipliers | Binary bilinear, mixed quadratic Boolean, and GF4-linear tensor searches were checked algebraically. | None found a replacement below 42 binary ANDs. The latest mixed-quadratic and GF4 searches each tried 64 million proposals; neither proves optimality. |

Evidence and exact scopes are collected in
[PHASE_AND_GROUPED_RESEARCH.md](PHASE_AND_GROUPED_RESEARCH.md),
[GROUPED_FACTOR_RESEARCH.md](GROUPED_FACTOR_RESEARCH.md),
[LOOKAHEAD_RESEARCH.md](LOOKAHEAD_RESEARCH.md),
[the rejected thin-update result](results/deferred_padding_thin_rejected.json),
and [the Goldschmidt checkpoint](results/gold_row_checkpoint.json).

**What the bounds establish.** For the fixed native chemistry shifts and gap
floor, the effective normalization must exceed 46.78617959252637. The selected
48.59417497035283 therefore has less than 3.87% improvement available from
normalization alone. A separate physical matrix-rank witness requires a
uniform factor-width limit of at least 10: at least one physical factor must
have matrix rank 10 or greater. This rules out the proposed B2/B4/B8
separate-factor families even after changing the factor-column gauge.
Reference-spectrum residual bounds also exclude total physical factor rank
275 or smaller at the 20 mHa threshold.

These bounds do **not** bound overall competition performance. They do not
rule out better lookup circuits, different factorizations, a different decoder,
or a different inverse algorithm. Their useful effect is to avoid spending more
time on the specific normalization-only and undersized-factor proposals.
See [NORMALIZATION_BOUND.md](NORMALIZATION_BOUND.md),
[normalization evidence](results/normalization_global_bound.json),
[matrix-rank evidence](results/physical_factor_rank_bound.json), and
[spectrum residual bounds](results/rank_error_lower_bounds.json).

**The remaining lead at shutdown.** The final experiment revisited the padding
cost used to reject a Wronskian-based decoder. Controlled affine transvections
reduce the N70 squarefree cofactor-padding adapter from 28,069 to **7,140 Toffolis
per direction**. The completed quantum test passes 4,096 cases each forward,
inverse, forward/inverse composition and fresh inverse/forward composition,
with arbitrary coefficients and arbitrary flags. It also checks phase and
clean scratch. The compiled component has 875 public qubits and no additional
live quantum workspace; this is not a complete decoder width.

That result invalidates the older padding-cost premise for rejecting this
route. It does not establish a submission improvement. A complete implementation
still needs a suitable inverse-only phase kernel, the padded-polynomial/H
adapter and a total reversible Sylvester transform. The existing 62,098T
Sylvester prototype was only checked on generic degree-17 polynomials and
explicitly excludes degree gaps. No total replacement was implemented before
the wrap-up request.

The tested padding implementation is [transvection_repad.py](transvection_repad.py),
with evidence in [transvection_repad_n70_check.json](results/transvection_repad_n70_check.json).
It is preserved outside the submission.

**Assessment.** The target remains a substantial distance away. From the current
score, reaching 14.0 requires another 3.338-fold reduction in the product of the
three task costs. If spread evenly, that means about a 1.494-fold cost reduction
in each task. The observed gains do not support a claim that small refinements
will close that gap.

The experiments support pursuing a larger representation or algorithm change
if work resumes. They do not establish that 14.0 is impossible. The competitor's
reported 16.6 result was not independently reproduced or available as an
implementation for comparison.

The best submission, native results, source manifests, reproducible builders,
rejected examples and research logs are retained. No host verifier or scoring
code was modified. [README.md](README.md) contains the selected rebuild and
verification commands; [CONTINUATION.md](CONTINUATION.md) records the detailed
handoff state.
