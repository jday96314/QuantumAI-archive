# Optimized circuit submission

The best submission found in the one-hour run has a **local overall score of 2.254539**, compared with 1.0 for the host baselines. The ready-to-upload archive is [submission.zip](submission.zip); its unpacked files are in [submission/](submission/).

All three unchanged host verifiers passed at their default sizes and shot counts, using an independent validation salt. The geometric mean uses the recorded local baseline costs. Official leaderboard salts and reference costs are private, so this is a local score.

| Benchmark | Baseline cost | Verified cost | Improvement | Scored qubits |
| --- | ---: | ---: | ---: | ---: |
| `invsqrt` | 2,076,774.233715 | 329,225.286361 | 6.3081× | 123 |
| `dqi` | 24,439,463,583,672.929688 | 15,548,775,312,944.798828 | 1.5718× | 1,861 |
| `dfthcsos` | 3,355,599,588.315537 | 2,903,278,577.426011 | 1.1558× | 739 |

## Inverse square root (`invsqrt`)

The final cost is **329,225.286**, down **84.15%** from 2,076,774.234: a **6.3081× improvement**. Width fell from 260 to 123 qubits, and raw Toffolis fell from 7,735 to 2,426. The scored cost also includes rotation synthesis.

### What helped

The largest gain came from replacing the baseline's cubic initial estimate and Newton refinement with a **direct piecewise quadratic** over 384 normalized mantissa intervals. This removed the Newton square and two large multiplications. Each interval supplies three integer coefficients through QROM; Horner arithmetic evaluates the quadratic using only the precision needed at each stage.

Subsequent changes reduced both gates and live workspace:

- Binary normalization and denormalization shifts replaced the ten unary shift flags with four binary flags.
- Measurement-based lookup deletion cleared the coefficient table without repeating its full forward lookup. Clean high bits of the output register also served as temporary workspace.
- Registers were trimmed to the coefficients' actual bit lengths. A custom unequal-width adder used temporary ANDs and measurement to erase carries, allowing shorter partial products without allocating zero padding.
- Exhaustive interval searches adjusted nearby integer slope and curvature coefficients, then fitted the constant term to the actual truncated arithmetic. The final fractional precisions are **28, 19, and 8 bits**. Arithmetic was checked over every normalized input before full circuit verification.

The following are cumulative, fully simulated candidates. Multiple changes sometimes landed together, so these measurements do not isolate every change's individual contribution.

| Candidate | Raw task cost | Qubits |
| --- | ---: | ---: |
| Host baseline | 2,076,774 | 260 |
| Table-assisted linear estimate plus Newton | 1,452,949 | 223 |
| First direct quadratic | 592,228 | 167 |
| Precision search, borrowed workspace, lookup deletion | 480,792 | 149 |
| Binary shifts | 425,681 | 143 |
| Trimmed coefficient registers | 403,372 | 141 |
| Custom short adder | 380,485 | 133 |
| Final integer-fitted quadratic | **329,225** | **123** |

Sources: [linear estimate](results/linear-default.json), [first quadratic](results/poly-default.json), [precision search](results/poly-search-0.json), [binary shifts](results/poly-binary.json), [register trimming](results/poly-trim.json), [short adder](results/poly-short.json), and [final verification](results/final-validation/invsqrt.json).

### What did not work

- The linear-estimate/Newton design passed, but retaining Newton's expensive arithmetic limited its improvement. Cubic and quartic interpolation candidates in the precision searches also cost more than the selected quadratic.
- Dynamically allocating shorter partial products while using kickmix's built-in unequal-width adder reduced width to 134 qubits but increased Toffolis to 3,180. Its [verified cost was 458,922](results/poly-dynamic.json), worse than the 403,372 candidate it was intended to improve. This motivated the custom adder.
- The tested 27-bit output configurations exceeded the relative-error tolerance. Further reducing the final curvature precision to seven bits reached **7.615e-8** worst-case relative error, above **2^-24 ≈ 5.960e-8**, and was rejected before circuit verification. See the [integer fitting log](results/integer-polynomial-search.log) and [last precision attempt](results/integer-polynomial-last-search.log).

## Reed–Solomon decoder (`dqi`)

The final cost is **15.549 trillion**, down **36.38%** from 24.439 trillion: a **1.5718× improvement**. Mean Toffolis fell from about 5.725 million to 3.769 million; width fell from 1,893 to 1,861 qubits. Both matter because this task scores `repetition_factor * Q^2 * mean_Toffolis`.

### What helped

The decoder retains the host's **222 synchronized EEA steps** and its overall decoding algorithm. The main changes are exact finite-field arithmetic replacements:

- Binary polynomial CRT first reduced multiplication from the host's 63 Toffolis to 47, then 45, using smaller field products and free binary linear maps. Extending this substitution into division and inversion chains improved the whole decoder.
- A quadratic subfield tower, **GF(4096) → GF(64) → GF(8)**, replaced the long inversion chain. The norm identity reduces a large inverse to smaller multiplications and one smaller inverse. GF(8) inversion is a quadratic Boolean map requiring only three ANDs.
- Evaluation and interpolation over **GF(4)** reduced a GF(4096) multiplication to **42 Toffolis** and a GF(64) multiplication to 15. Together with the tower, the final GF(4096) inverse costs **66 Toffolis**, versus the host's 315, while retaining only 18 scaffold qubits.
- Seven-bit metadata counters saved eight more qubits in the complete circuit and slightly reduced arithmetic cost.

These cumulative candidates all retained 222 iterations and passed 102,400-shot verification with zero errors:

| Candidate | Raw task cost, trillions | Qubits |
| --- | ---: | ---: |
| Host baseline | 24.439 | 1,893 |
| Binary CRT in the multiplication sweep | 22.226 | 1,893 |
| Binary CRT throughout arithmetic chains | 20.030 | 1,896 |
| Subfield tower inversion | 16.354 | 1,869 |
| GF(4) multiplication construction | 15.752 | 1,869 |
| Compact metadata counters | **15.549** | **1,861** |

Sources: [sweep substitution](results/dqi-crt-222.log), [full arithmetic substitution](results/dqi-full-crt-222.log), [tower inversion](results/dqi-tower-222.log), [GF(4) construction](results/dqi-fast-222.log), and [compact counters](results/dqi-compact-222.log). The final artifact passed a further [independent-salt run](results/final-validation/dqi.json). Arithmetic checks covered all bilinear basis pairs, all 4,096 inverse inputs and divisors, restored workspace, and phase corrections; see [arithmetic test results](results/final-arithmetic-tests.log).

### What did not work

Reducing the EEA iteration budget produced a poor cost/error tradeoff. [190 iterations](results/dqi-190.log) and [214 iterations](results/dqi-crt-214.log) had phase errors on roughly half the sampled shots and failed verification. Shorter runs close to the baseline were more deceptive: a native 218-step decoder had zero errors in an 8,192-shot check, but a later 218-step arithmetic variant had **17 phase errors in 102,400 shots**. A 216-step variant had **32**. The variants differ in arithmetic and sampling seeds, so these runs do not establish an exact failure rate for the iteration count alone.

Although the 216- and 218-step variants passed the verifier's allowed-error gate, their repetition penalties increased to about **1.2804** and **1.2476**, compared with **1.1913** for zero observed errors at the default shot count. Their costs of 21.078 and 21.200 trillion were worse than the full-length arithmetic candidates. The final submission therefore keeps all 222 iterations. See [the 216-step run](results/dqi-full-crt-216.log) and [the 218-step run](results/dqi-full-crt-218.log).

## Chemistry block encoding (`dfthcsos`)

The final cost is **2.903 billion**, down **13.48%** from 3.356 billion: a **1.1558× improvement**. This task scores `lambda_eff * T_total * Q_total`, so the work combined circuit changes with classical factorization changes. The final scored width is 739 qubits, slightly above the baseline's 730; lower gate cost and normalization more than offset that increase.

### What helped

The **128-qubit ANGLE-circuit cap** constrained lookup batching. The final QROM loads nine 12-bit words per batch. It computes high address prefixes separately and replaces the last unary-iteration ancilla with a multi-target CCX: CNOT conjugation lets one Toffoli toggle an entire output vector. This fits nine batches into exactly 128 qubits per circuit. Measurement-based prefix deletion then reduced the executed lookup cost further.

For the classical factorization, orthogonal rotations of ensemble columns preserve the represented tensors while changing the LCU one-norm. Pairwise coordinate searches reduced that normalization. A matrix-norm-weighted truncation then compressed **C=60 to C=29**, retaining **R=65, B=60**, and the original 12-bit angle words. The unrefitted 29-column candidate was just outside the accuracy limit. Refitting its nonnegative squared column strengths brought it inside; identity offsets were scaled inversely to preserve the one-body correction. A final orthogonal refinement reduced the normalization again.

| Cumulative candidate | Raw task cost, billions |
| --- | ---: |
| Host baseline | 3.356 |
| Nine-word QROM plus initial ensemble rotations | 3.120 |
| Finer ensemble rotations | 3.105 |
| 36-column compression and rotations | 2.986 |
| Measurement-based prefix deletion | 2.945 |
| Refitted 29-column factorization and final rotations | **2.903** |

Sources: [nine-word lookup](results/chemistry-j9.json), [rotation refinement](results/chemistry-refined.json), [36-column factorization](results/chemistry-c36.json), [prefix deletion](results/chemistry-j9-eq.json), and [final verification](results/final-validation/dfthcsos.json).

### What did not work

- An ordinary nine-word QROM would require 131 qubits, exceeding the 128-qubit cap. The pebbled construction was needed to fit nine words.
- Unweighted SVD truncation did not reliably preserve accuracy or improve normalization. At **C=30**, its two-body error was **22.237 mHa**. Even the passing C=40 tensor candidate raised `lambda_eff` to about **51.38**, versus the baseline's 49.30. This motivated weighting compression by the decoded matrix norm and then rotating the retained columns. See [compression measurements](results/chemistry-compression.log).
- Reducing the orbital-vector basis B using weighted pivoted QR and least squares was much too destructive: even **B=58** produced **118.47 mHa** two-body error, far above the 20 mHa limit. See [basis compression measurements](results/chemistry-basis-compression.log).
- Weighted **C=28** truncation had about **20.56 mHa** two-body error. Strength refitting lowered it only to **20.316 mHa**, so it remained invalid. The unrefitted C=29 candidate also failed at **20.0186 mHa**; refitting C=29 succeeded at **19.3953 mHa**. See [the C=28 refit](results/chemistry-refit.log), [the C=29 rejection](results/chemistry-c29-weighted.json), and [the successful refit](results/chemistry-refit29.log).
- A further direct effective-normalization optimization terminated abnormally and returned the same metrics, so it supplied no improvement. See [that optimizer log](results/chemistry-effective29.log).

The final chemistry tensor errors are 3.257691 mHa and 19.395341 mHa, below the 20 mHa limits. Lambda is 292.692404468 and E_gap is 3.980013438, above their strict acceptance floors.

## Reproduction and validation

The canonical builders are [build_submission.py](build_submission.py), [polynomial_invsqrt.py](polynomial_invsqrt.py), [short_adder.py](short_adder.py), [fast_dqi.py](fast_dqi.py), [gf4_multiplier.py](gf4_multiplier.py), [tower_inverse.py](tower_inverse.py), [compact_dqi.py](compact_dqi.py), and [chemistry_qrom.py](chemistry_qrom.py). Exact fitted parameters are frozen in [parameters/](parameters/); optimizer reruns are unnecessary to rebuild the circuits.

Run from the repository root with Python 3.13 and a compatible C++ runtime. The provided environment script points to the environment used here (NumPy 2.5.4 and SciPy 1.18.1); update its interpreter/runtime paths on another machine.

```bash
source initial_experiments/env.sh
"$EXP_PYTHON" initial_experiments/build_submission.py --output /tmp/quantumai-optimized-rebuild
"$EXP_PYTHON" initial_experiments/verify_submission.py --submission /tmp/quantumai-optimized-rebuild --output /tmp/quantumai-optimized-validation --salt independent_reproduction_validation_20261023 --workers 4
```

The authoritative measurements are [summary.json](results/final-validation/summary.json) and the three adjacent verifier outputs. [metric-contract.json](results/metric-contract.json) records checks with the host metric's actual circuit and NPZ readers, including instruction allowlists, widths, shapes, and regular files. [archive.json](results/archive.json) records the ZIP size, CRC validation, and SHA-256. The exploratory scripts and logs retain unsuccessful precision, iteration, and factorization experiments.

Archive SHA-256: `3a413fcaef90d4d39055d4c305a1978be95cdd8b7eb51e53d6e46edc0f9287e1`.

Host files and the pre-existing root README changes were preserved.
