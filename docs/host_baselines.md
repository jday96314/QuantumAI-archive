# How the host baselines work

The supplied [sample submission](../from_host/sample_submission) contains concrete reference circuits and a chemistry factorization. These are starting implementations, not lower bounds or requirements to retain their algorithms. The unchanged host verifiers define the accepted behavior; alternative algorithms can replace the references.

Recorded full-default local measurements are below. These costs normalize to a combined score of 1.0 against themselves. Official baseline denominators are measured under their own split salts, so these are local reference values.

| Task | Raw cost | Scored qubits | Mean Toffolis | Verified shots |
| --- | ---: | ---: | ---: | ---: |
| Inverse square root | 2,076,774.233715 | 260 | 7,735 | 200,540,928 |
| DQI | 24,439,463,583,672.93 | 1,893 | 5,725,108.3821 | 102,400 |
| Chemistry | 3,355,599,588.315537 | 730 | 93,237 | 254,464 |

See [baseline profiling](baseline_profiling.md) for measurement settings, stage attribution, and runtime. [Verifier details](verifiers.md) explains acceptance and score formulas.

## Inverse square root: a polynomial seed followed by Newton refinement

The [reference builder](../from_host/verifiers/invsqrt/nr_reference_circuit.py) applies a phase proportional to `1/sqrt(x)` while restoring its 22-bit input and clearing all temporary values.

1. Count leading zero pairs in a unary register and shift x by powers of four. The normalized mantissa M lies in `[1,4)`, with the shift recording the scale needed afterward.
2. Choose one of four cubic polynomials over mantissa intervals `[1,1.5)`, `[1.5,2)`, `[2,3)`, and `[3,4)`. Coefficients are loaded by a small Clifford/AND network. Truncated fixed-point Horner arithmetic produces an initial inverse-square-root estimate y0.
3. Refine with the Newton expression `p = y0 × (3 − M×y0²)`, with specified partial-product truncation. This approximates `2/sqrt(M)`; output scaling supplies the factor needed for `1/sqrt(x)`.
4. Apply bitwise phase rotations from the result, then reverse the arithmetic and normalization.

The fixed-point model mirrors truncation and checks overflow; substituting floating-point Newton arithmetic would not reproduce the circuit. Default fractional precisions are 14 bits for the interval offset, 17 for coefficients, 28 for the square and correction product, and 29 for the final refined product.

Arithmetic dominates: Horner costs 2,142 Toffolis, multiplying M by the squared estimate 2,086, the final correction product 1,734, and the initial square 884, including reverse work. Together these consume 88.51% of raw Toffolis. Rotation synthesis raises 7,735 raw Toffolis to 7,987.5932 equivalents, which multiplied by 260 gives the score. The direct quadratic replacements in [iteration 1](optimization_iteration_1.md) avoid the Newton square and products.

## DQI: synchronized reversible Euclidean decoding

The [reference builder](../from_host/verifiers/dqi/dqi_reference_circuit.py) implements the Proos–Zalka-style synchronized extended Euclidean algorithm over GF(4096). From 70 syndrome words it computes error-locator and error-evaluator polynomials, scans all 4,095 nonzero field positions, and reverses the decoder.

The forward pass uses **222 fixed synchronized iterations** for the benchmark. Metadata, controlled shifts/swaps, and multiply-add sweeps implement degree-dependent polynomial updates reversibly. The Chien search finds roots of the locator; Forney evaluation recovers error values and contracts them with classical measurement masks into the required sign. The decoded error vector is not retained as output. Reverse EEA restores the input and removes polynomial workspace.

The circuit contains 840 syndrome qubits and 1,053 scratch qubits. Forward and reverse EEA each cost about 2.109 million Toffolis in the profiling batch, together 73.69% of the total. Chien/Forney costs about 1.444 million, or 25.23%; alignment contributes the remainder. Polynomial position sweeps are the largest arithmetic target.

The full recorded run has zero value and phase errors, but its statistical upper confidence bound still gives a repetition factor of **1.19126046**. Its cost is that factor times `1,893²` times the executed mean Toffolis. Static gate counts differ slightly because classically conditioned execution varies by shot. Later experiments first replaced field arithmetic, then replaced the synchronized decoder itself.

## Chemistry: an existing factorization plus ordinary angle QROMs

The sample uses a **65-ensemble, 60-vector, 60-column** factorization: `(R,B,C) = (65,60,60)`. Its log-depth Givens tree represents each 76-dimensional unit vector using 75 angles. There are `76 + 65×60 = 3,976` rows: 76 one-body rows followed by the factorization vectors.

Angles are quantized to **12 bits**. Ten ordinary unary-iteration QROM circuits load eight words at a time, padding the final batch. Each oracle has 119 qubits and 3,974 Toffolis. The [example builder](../from_host/verifiers/dfthcsos/example_submission/angle_k_reference_circuit.py) regenerates these gates from the sample NPZ's extra `angle_words` array; that convenience array is not the verifier's authoritative source of words.

The submitted lookup circuits implement the angle-loading portion, while the host analytically prices the rest of block encoding:

| Component | Toffolis |
| --- | ---: |
| Forward and reverse ANGLE lookups | 79,480 |
| Inner state preparation | 9,461 |
| Outer state preparation | 846 |
| Givens rotations | 3,450 |
| Total | **93,237** |

Width is **186 persistent + 544 shared temporary = 730**. The persistent part includes 96 angle-word qubits and 90 preparation qubits; inner preparation sets the temporary maximum. The ten individual oracle widths are not summed. This explains why reducing lookup ancillas alone need not lower the scored width, while batching can still reduce gates.

Decoded words produce one-body error about **3.1852 mHa** and two-body error **18.6396 mHa**. Raw normalization is **308.45594**, SOS gap **3.965477 Ha**, and effective normalization **49.301369**. The score multiplies that effective normalization by 93,237 and 730. Quantization changes the normalization from the unrounded factorization's value and leaves limited physical-error headroom.

The fixed symmetry shifts and checksum-validated FCIDUMP are supplied by the verifier. Optimization can change weights, identity offsets, dictionary shape, tree, angle precision, and oracle implementation, subject to the contract. [Iterations 1](optimization_iteration_1.md), [2](optimization_iteration_2.md), and [3](optimization_iteration_3.md) show why fewer factor columns, fewer ensembles, and a shared dictionary are distinct optimizations with different accuracy and normalization trade-offs.
