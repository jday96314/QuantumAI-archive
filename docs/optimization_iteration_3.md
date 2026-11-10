# Optimization iteration 3

The [third campaign](../initial_experiments/3/EXPERIMENT_REPORT.md), including its later controlled sign-mixing study, reached a local combined score of **9.0096229084**, up from iteration 2's **4.5267**. Its pre-follow-up checkpoint scored 9.0085406883. Some builders were inherited; the table compares selected outcomes rather than assigning every idea's origin to this directory.

| Task | Iteration 2 cost | Final cost | Improvement over host |
| --- | ---: | ---: | ---: |
| Inverse square root | 311,390 | 190,711.9555 | 10.8896× |
| DQI | 3.268 trillion | 1.612737419463 trillion | 15.1540× |
| Chemistry | 1.804 billion | 757,163,489.50 | 4.4318× |

## Inverse square root: centered interpolation and tighter reversible scheduling

Binary normalization with parity folded into lookup, centered interpolation, and integer fitting only on reachable mantissas reduced table/arithmetic cost. Short accumulators delayed carry integration; rounded partial products reduced coefficient widths. Radix-four Booth accumulation shared work between adjacent coordinate bits.

Hybrid carry schedules, borrowed output bits, and delayed allocation reduced live width. Measuring unused high curvature bits early was beneficial only when their phase correction was deferred until coefficient reconstruction during reversal. Private-wire compaction completed the reduction to **94 qubits**.

Verified cumulative costs fell from 298,067 for binary-centered interpolation to 233,306 with rounded products, 209,397 with adaptive Booth arithmetic, 194,025 at 98 qubits, and **190,712** after delayed phase repair. The final circuit had about 1,778 mean raw Toffolis, 34 rotations, and 2,028.851 Toffoli equivalents. It passed **200,540,928 shots**.

Dropping a low rotation passed but worsened remaining synthesis costs. Alternate rounding/tie rules and global phase offsets did not improve the tested grids. Immediate curvature phase repair added too much lookup work; Boolean synthesis of high curvature bits was also uncompetitive. Wider low-gate circuits could lose once width and synthesis margin were priced.

## DQI: consume syndrome storage rather than retaining every register

An attributed adaptation of the public `andigu/uncompute` Berlekamp–Massey decoder supplied a different reversible framework, runtime masks, and host-field basis conversion. Consumed syndrome registers became storage for the backup polynomial bank. Parked inputs and reversible stop guards handled exceptional leading coefficients.

Recursive subfield phase operations produced an **80-Toffoli** point kernel. Selecting a six-bit inverse norm before the final product avoided full-width intermediate registers; using an armed locator coefficient directly as a scale removed another duplicate.

| Verified stage | Cost, trillions |
| --- | ---: |
| Initial public-decoder adaptation | 1.829622 |
| Syndrome-consuming stream3 | 1.772682 |
| Corrected final-round guards | 1.681616 |
| Compact selection and scale | **1.612737** |

The selected circuit used **1,429 qubits**, **645,755 Toffolis**, and repetition factor **1.223014**. Full verification recorded **zero value errors and eight phase errors in 102,400 shots**.

Early journal replay, translation, and iterator bugs broke reversal and were fixed. A stale stop guard rejected valid weight-35 paths. Shrinking the bank without matching consumption to its lifetime produced roughly 1,400 phase errors; extra discrepancy gating added about 32,000 Toffolis without fixing the cause. Some compaction schedules created allocation conflicts. Small-instance passes were diagnostics only; selection used all 4,095 evaluation points.

## Chemistry: a shared dictionary and direct integer-angle fitting

Early compression of the inherited 35-ensemble shape reached about 1.758 billion. A global shared dictionary then reduced angle-table rows substantially. PSD projection had to respect the symmetric feature Gram metric: naive Euclidean clipping could turn a roughly 19 mHa fit into 68 mHa.

Integer-angle moves, rank-constrained weight refits, and dictionary growth/reseeding finally produced a valid **eight-bit, 719-vector, 350-column** dictionary. Column alignment/rotations, identity-offset fitting, and constrained whitened-coordinate updates reduced normalization while preserving or carefully spending the physical-error budget. One-body angle words were separately optimized and regenerated after weight changes.

| Verified stage | Cost, millions |
| --- | ---: |
| Polished shared dictionary, 11 bits | 909.956 |
| 10-bit fit | 838.859 |
| Initial valid eight-bit dictionary | 823.241 |
| Constrained refinement | 757.436 |
| Weak mixing-penalty follow-up | **757.163** |

The selected preparation used **29,440 Toffolis and 526 scored qubits**. Six thirteen-word lookup batches fit under the individual 128-qubit cap. More aggressive compression, nonnegative weights, straight-through angle fitting, and grouped blocks did not produce better valid submissions. One grouped float fit reached 14.068 mHa but about 31 mHa after eight-bit quantization. Published THC factors helped initialization but did not directly meet the shifted reference accuracy. A wrong spectral ordering also invalidated an accelerated fit and was excluded.

## Controlled follow-up: sign mixing was a weak surrogate

The study tested `M(w) = Σc [||w_c||₁² − (Σb w_bc)²]` with fixed dictionary, tree, angle precision, and two-body words. Matched 16,000-step weight/offset refits separated penalty effects from simply continuing optimization.

The previous normalization was 48.912811; an unpenalized control reached **48.895752**, while the selected weak error-scaled penalty reached **48.895187**. Chemistry improved **0.036031%** over the old checkpoint, but only **0.001154%** beyond the matched control. Stronger penalties lowered mixing while increasing the gap enough to worsen effective normalization and cost. Matched physical-column rotation controls reached the same conclusion: less cancellation alone did not guarantee a lower objective.

Seven new chemistry candidates passed 50,880-shot native verification. The selected physical errors were **11.9968 mHa one-body** and **19.9953 mHa two-body**, with gap **3.915138 Ha**, close to the acceptance floor. The archive audit bound selected files to native results; unchanged inverse/DQI reused existing full verification. [The detailed report](../initial_experiments/3/EXPERIMENT_REPORT.md) contains controls, coefficient conventions, candidate links, and provenance.

The main lesson was to evaluate actual quantized physics and the complete objective, and to credit continued fitting separately from a proposed regularizer. [Iteration 4](optimization_iteration_4.md) tested further storage, lookup-representation, and decoder changes.
