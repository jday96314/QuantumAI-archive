# Optimization iteration 2

The [second campaign](../initial_experiments/2/README.md) raised the local combined score from **2.2545 to 4.5267**. Its main gains came from a different DQI decoder and fewer chemistry ensembles. Selected measurements use unchanged host verifiers at default sizes and shots.

| Task | Iteration 1 cost | Selected cost | Improvement over host |
| --- | ---: | ---: | ---: |
| Inverse square root | 329,225 | 311,390 | 6.67× |
| DQI | 15.549 trillion | 3.268 trillion | 7.48× |
| Chemistry | 2.903 billion | 1.804 billion | 1.86× |

## DQI: generic fixed-position Euclid

The key change replaced 222 synchronized EEA iterations with up to **35 generic degree-one quotient updates**. Two 71-slot arrays combine remainder and cofactor storage. Slots vacated by zeroed remainder coefficients become new cofactor slots, with positions known at circuit construction time. This removes degree-register routing, full shifts, and array swaps. Quotients are erased using the cofactor relation, and an active-step counter handles early completion.

Degree-one quotients are a typical-case assumption, not a universal property. Accidental vanishing leading coefficients produce rare exceptional cases; the final run reported **359 phase errors and zero value errors in 102,400 shots**, increasing the repetition factor to about **1.68**. The gate reduction outweighed that penalty.

| Version | Mean Toffolis | Qubits | Cost, trillions |
| --- | ---: | ---: | ---: |
| Previous synchronized decoder | 3,768,745 | 1,861 | 15.549 |
| Cheaper Forney kickback | 3,715,522 | 1,861 | 15.329 |
| Generic Euclid | 819,571 | 1,832 | 4.643 |
| Evaluation/interpolation updates | 616,710 | 1,832 | 3.490 |
| In-place sums and counter flags | 627,702 | 1,780 | 3.333 |
| Measurement-erased quotients and cleanup scheduling | about 621,900 | **1,766** | **3.268** |

Evaluation/interpolation made a linear-polynomial update use roughly 73 field products per step instead of 142; constant field maps are CNOT networks. The search kernel used 89 Toffolis per point: 66 for inversion, 12 for mask parity, and 11 for a zero test. Measurement cleanup avoided explicitly reversing all product work.

Repairing exceptional degree gaps was identified as a possible future gain, but not built in this iteration. The passing result therefore remains a priced approximate decoder. Fewer Toffolis alone did not select the winner: the lean variant deliberately spent gates to reduce squared width cost.

## Chemistry: refit fewer ensembles, then quantize and rotate

GPU fitting optimized the physical two-body tensor, with TF32 disabled, using full per-ensemble column capacity **C=B**. This differs from iteration 1's compression of C while keeping 65 ensembles. Vectors were encoded as 12-bit Givens words, weights refitted with those vectors fixed, offsets optimized, and one-body eigenvectors encoded separately. Orthogonal rotations of weights and offsets preserved physics while lowering normalization.

| Shape R/B/C | Quantized two-body error | Cost, billions |
| --- | ---: | ---: |
| 55/60/60, rotated | 14.92 mHa | 2.654 |
| 50/60/60, rotated | 15.02 mHa | 2.431 |
| 40/60/60, rotated | 18.50 mHa | 2.018 |
| **35/60/60, rotated** | **19.60 mHa** | **1.804** |

The selected fit used 49,404 Toffolis and effective normalization about 49.29. The earlier nine-word lookup construction was reused.

Aggressive shapes failed after quantization: R30/B60/C60 reached 20.2 mHa, R35/B55/C55 21.1 mHa, and tested B≤50 configurations stalled at ≥56 mHa. Reducing the tested fit to 11-bit angles gave 51.7 mHa. This was a failure of those fits/quantization choices, not a universal eleven-bit impossibility. Re-optimizing offsets after rotation produced no useful selected improvement.

## Inverse square root: cheaper shifts and flag deletion

The quadratic design was retained. Zero-aware controlled shifts replaced cyclic rotations, saving 126 Toffolis. Measurement-based deletion of normalization flags with conditional phase-if-zero repair saved about another 11 executed Toffolis on average. Cost fell from 329,225 to **311,390**.

## Evidence and handoff

The [iteration README](../initial_experiments/2/README.md) names the selected verifier outputs, builders, fit scripts, and package assembly. The main outcomes were a large typical-case decoder gain with explicitly priced rare failures, and a smaller physical factorization found through refitting rather than truncation alone. [Iteration 3](optimization_iteration_3.md) moved to syndrome-consuming decoding and a shared chemistry dictionary.
