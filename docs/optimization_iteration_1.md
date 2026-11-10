# Optimization iteration 1

This is the campaign in [initial_experiments/1](../initial_experiments/1/README.md). It improved the local combined score from **1.0 to 2.254539** in the recorded one-hour run. All selected artifacts passed unchanged host verifiers at default sizes/shots and an independent validation salt. Figures are cumulative candidate results, not isolated causal measurements.

| Task | Host cost | Selected cost | Improvement over host |
| --- | ---: | ---: | ---: |
| Inverse square root | 2,076,774 | 329,225 | 6.3081× |
| DQI | 24.439 trillion | 15.549 trillion | 1.5718× |
| Chemistry | 3.356 billion | 2.903 billion | 1.1558× |

## Inverse square root

The largest gain replaced the host's cubic seed plus Newton step with a **direct piecewise quadratic over 384 normalized intervals**. It removed the square and two large Newton products. Coefficients loaded through QROM feed precision-trimmed Horner arithmetic.

Binary shifts replaced unary normalization flags. Measurement-based table deletion avoided replaying the full lookup; idle high output bits supplied scratch. Unequal-width temporary-AND adders eliminated zero padding, and exhaustive integer fitting optimized coefficients against the actual truncated products. Selected fractional precisions were 28/19/8 bits. Width fell from 260 to **123**, and raw Toffolis from 7,735 to **2,426**.

| Cumulative stage | Raw cost |
| --- | ---: |
| Linear estimate plus Newton | 1,452,949 |
| First direct quadratic | 592,228 |
| Precision search and lookup/workspace changes | 480,792 |
| Binary shifts | 425,681 |
| Custom short adder | 380,485 |
| Integer-fitted quadratic | **329,225** |

Linear/Newton and cubic/quartic alternatives were more costly. A built-in unequal-width-adder variant saved width but increased gates, scoring 458,922 versus the preceding 403,372. Tested 27-bit outputs violated accuracy; seven-bit curvature reached relative error 7.615e-8 versus the 2^-24 limit of about 5.960e-8. Those precision reductions were excluded.

## DQI

The selected decoder retained all **222 synchronized host EEA iterations**, replacing exact finite-field primitives. Binary polynomial CRT reduced multiplication from 63 to 47 then 45 Toffolis. A GF(4096)→GF(64)→GF(8) tower made inversion cheaper; GF(4) evaluation/interpolation then reduced multiplication to **42 Toffolis**. The final GF(4096) inverse cost **66**, versus 315 in the host arithmetic. Seven-bit metadata counters saved additional width.

Cumulative costs fell through 22.226, 20.030, 16.354, and 15.752 trillion to **15.549 trillion**. Final width was **1,861** and mean Toffolis about **3.769 million**. The full 102,400-shot selected runs had zero errors; primitive checks covered basis-pair products, every inverse/divisor, phase, and clean scratch.

Shortening EEA was unsuccessful: 190/214-step tests had about half their phases wrong. A short 218-step diagnostic saw no errors, but longer tests of arithmetic variants found 17 errors at 218 steps and 32 at 216. Some passed acceptance but paid repetition factors around 1.2476/1.2804 and cost more. Different variants and seeds prevent attributing exact failure rates to iteration count alone.

## Chemistry

Nine-word batches of 12-bit angles initially exceeded the 128-qubit oracle cap. A pebbled prefix lookup and multi-target CCX construction fit them at exactly 128 qubits; measurement-based prefix cleanup saved further gates. Orthogonal column rotations reduced normalization without changing the represented tensors.

Matrix-norm-weighted compression reduced columns from 60 to **29**, keeping **R65/B60** and the original angle precision. Refitting nonnegative squared column strengths rescued the almost-valid truncation; inversely scaling identity offsets preserved the one-body correction. Final one-/two-body errors were **3.2577/19.3953 mHa**, with effective normalization about 48.10. Total scored width increased to **739**, but gate and normalization savings won overall.

Unweighted SVD was a poor physical proxy: C30 gave 22.237 mHa, while a passing C40 fit increased normalization. Reducing B to 58 gave 118.47 mHa. Weighted C28 still failed after refitting at 20.316 mHa; C29 improved from 20.0186 to 19.3953 and passed. A later normalization optimizer terminated abnormally and supplied no improvement.

## Evidence and handoff

The [detailed report](../initial_experiments/1/README.md) links cumulative candidates, failed fits, arithmetic tests, final default verifications, and archive/metric-contract checks. Frozen parameters and builders preserve the selected output. The main lesson was to replace expensive algorithms where possible and evaluate precision, lookup, and workspace together. [Iteration 2](optimization_iteration_2.md) subsequently replaced the decoder architecture and refitted the chemistry shape.
