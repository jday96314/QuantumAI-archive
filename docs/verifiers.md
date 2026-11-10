# Verifier and scoring details

This supplements the [competition overview](competition_overview.md) with details checked against the bundled [metric](../from_host/metric.py) and verifier source. It describes this repository's host snapshot, not a live leaderboard configuration. The [baseline guide](host_baselines.md) explains the supplied algorithms; [baseline profiling](baseline_profiling.md) records measured costs and evaluator runtime.

## Shared accounting and reproducibility

Raw task costs are lower-is-better. The combined score is `(baseline_inv/cost_inv × baseline_dqi/cost_dqi × baseline_chem/cost_chem)^(1/3)`, so higher is better overall. Official reference costs and salts are split-specific; local ratios against the recorded host runs are local scores.

Toffoli-like gates include CCX, CCZ, and CSWAP. The simulator counts executed gates, including classical conditions, so the scored mean can be below a circuit's static operation count. Clifford operations and classical control can still dominate simulation time even when they are cheap in the quantum objective. Runtime is not the leaderboard cost.

Seeds derive from canonical circuit/submission bytes combined with a salt. An explicit local seed overrides that derivation. Changing a circuit can change the sampled inputs, measurement outcomes, error count, and final score. Repeating identical bytes with the same settings and salt repeats the same samples; it supplies no independent evidence. Comments or archive repackaging do not necessarily change the canonical hash. Final comparisons should retain hashes, salts/settings, native verifier JSON, and the exact archive members.

## Inverse square root

The [verifier](../from_host/verifiers/invsqrt/verify_submission_invsqrt.py) expects one little-endian, 22-qubit register named `x`. Every other qubit is scratch. It checks restoration of `x`, zero scratch, and phase `exp(2πi/√x)` modulo one turn. Default inputs are **every integer from 1 through 3,133,452**, each repeated 64 times to test measurement branches: **200,540,928 shots**. Zero is outside the domain.

The tolerance at each input is `2^-24 / sqrt(x)` turns. Reference phases and tolerances use integer arithmetic and 64-bit fixed-point phase comparisons. Any observed restoration, scratch, or phase violation rejects the circuit.

The cost is `Q × (mean_Toffolis + total_T/4)`, with `Q = circuit.num_qubits`. The verifier charges explicit T angles and synthesis for non-Clifford rotations. For each synthesized rotation, `T(ε) = 0.56 log2(1/ε) + 5.3`; the shared allowed operator-norm error is chosen so the worst-case accumulated synthesis phase error fits the **smallest absolute margin** left by the simulated approximation. Consequently, deleting a rotation or reducing arithmetic precision can make the remaining rotations more expensive. A numerically passing approximation with no synthesis margin can still fail.

Full-domain coverage does not exhaust every possible measurement record. Shorter domains and fewer repeats are diagnostics, not interchangeable final scores.

## DQI

The [verifier](../from_host/verifiers/dqi/verify_submission_dqi.py) uses RS(4095,4025) over GF(4096): 70 twelve-bit quantum syndrome registers `S0`…`S69`, plus 4,095 twelve-bit **classical** outcome registers `E0`…`E4094`. The quantum input occupies 840 qubits; the outcome registers are not extra quantum input qubits. Registers must follow the exact schema and be disjoint.

Given an error vector `e` consistent with the syndromes and classical masks `M`, the circuit must restore the syndrome and clean scratch while applying sign `(-1)^Σ popcount(e_j & M_j)`. It need not return a classical decoded vector. Arbitrary `Z_POW` rotations are prohibited.

The default **102,400 shots** run in batches of 256. Error weights from zero through 35 are drawn from the squared principal-eigenvector weights of the verifier's DQI matrix, rather than uniformly. Positions, nonzero field values, and masks are sampled afterward. A decoder that works on generic degree-one Euclidean quotients can therefore have rare degree-gap failures that matter to scoring.

`value_errors` counts unrestored syndrome or dirty scratch; `phase_errors` counts wrong signs only on shots without value errors. Their sum avoids double-counting. Unlike the other tasks, a nonzero observed error count can pass: the one-sided Clopper–Pearson upper bound `ε+`, at confidence `1 − 2^-128`, must be **strictly below 1/16**.

The cost is `R(ε+) × Q² × mean_Toffolis`, where `R(ε) = ln(2) / [-ln(1/2 + 2√ε)]`. Even zero observed failures has a positive confidence bound: at default shots the recorded host run has `ε+ ≈ 0.00086606` and `R ≈ 1.19126`. Reducing shots raises this penalty even if no errors are seen. `PASS` means the priced error criterion was met, not exact decoding. Width savings must be evaluated together with gates and repetition; a cheaper kernel can score worse on its new sampled error count.

## Chemistry

The [scoring guide](../from_host/verifiers/dfthcsos/SCORING_README.md) and [verifier](../from_host/verifiers/dfthcsos/verify.py) define a hybrid physical/circuit contract. The system is FeMoCo76o113e: 76 spatial orbitals and 113 electrons.

The NPZ supplies `U (R,B,76)`, `W (R,B,C)`, `H2_I (R,C)`, scalar `H1_I`, `Hsym_lower (2926,)`, and `tree_edges (75,2)`. **Only U's shape supplies the orbital-vector layout; the physical vectors are decoded from circuit-loaded angle words.** Submitted `H1_I` and `Hsym_lower` are required by the file contract but ignored physically: the verifier uses its checksum-validated fixed symmetry shifts. Altering them does not optimize the Hamiltonian shift.

There are `76 + R×B` addresses and 75 Givens words per address. An integer word means `θ/π = word/2^(brot−1)`. Register names/layout encode selection, batch size J, and word precision. The ordered `angle_k_0`, … family must have consistent public registers. Each individual oracle is capped at **128 qubits**; this does not cap the scored full preparation width.

The family streams batches through the same word registers using **XOR differences**: oracle zero XORs in batch zero, and oracle k XORs `batch[k] XOR batch[k−1]`. Thus it replaces the preceding batch along the forward trajectory, and reversing the family clears the words. Loading each full batch independently would implement the wrong interface. `selection` has exactly `ceil(log2(76 + R×B))` qubits; extra address bits cannot hide scratch.

Every address receives 64 measurement samples. Circuits run forward and in reverse; the verifier recovers words and checks stable outputs, restored selection, zero workspace, no unwanted phase, zero padding, and reverse clearing. `Z_POW` is prohibited. Any observed circuit violation rejects the submission. These trajectories do not constitute an exhaustive arbitrary-incoming-target oracle test; supplemental random-target tests can catch mistakes beyond the native trajectory.

The recovered quantized words, tree, and weights determine physical residuals and normalization. Acceptance requires:

- One-body and two-body reconstruction residuals each **≤ 0.020 Ha**. A 12 mHa fitting target in an experiment is a local margin, not the host's one-body threshold.
- Raw `Lambda > 49.0` and `E_gap > 3.914898 Ha`, both strict.
- `0 < E_gap < 2 Lambda`, where `E_gap = E_upper − E_SOS` and fixed `E_upper = −22140.27335 Ha`.

The cost is `sqrt(E_gap × (2 Lambda − E_gap)) × T_total × Q_total`. There is no independent effective-normalization floor imposed by the verifier. The two reconstruction tests concern coefficient residuals; they do not certify a 20 mHa many-electron energy error.

`T_total` adds the simulated forward/reverse ANGLE gates to analytical outer/inner state preparation and Givens costs. `Q_total` adds persistent state-preparation and angle-word storage to the **maximum**, rather than sum, of shared temporary workspaces. Reducing an oracle's workspace alone helps width only if that workspace controls the maximum. Shape and precision changes can cross discrete QROAM thresholds and sharply increase scored width.

The FCIDUMP checksum, float32 rounding of two-electron integrals, and float32 symmetry-shift semantics are part of the reference calculation. Float fits, proxy losses, and unquantized arrays cannot replace native evaluation. Chemistry sampling does not report a pooled Clopper–Pearson certificate; 64 repeated branches per address are not a universal end-to-end error bound.

## Packaging checks and local use

The metric validates files before running task verification. A standalone verifier pass does not establish that an archive satisfies the metric's additional checks:

- Exactly one inverse circuit, one DQI circuit, and one `dfthcsos` directory; chemistry oracle indices start at zero without gaps or duplicate text/binary versions.
- Inputs must be regular files inside the extracted submission, without links or special files. All registers must be named; instruction and register-count limits apply. Qubit/classical-bit IDs must stay below `2**30 // 1536`.
- Required NPZ arrays must have consistent positive named dimensions and permitted plain numeric types, with no object/structured or big-endian arrays. The reader requires NPY 1.0 headers and caps the required arrays at **64 MiB** total.

Use the [root setup instructions](../README.md) for the bundled Python/native runtime, then invoke the unchanged verifiers:

```bash
python from_host/verifiers/invsqrt/verify_submission_invsqrt.py PATH/invsqrt.kmb
python from_host/verifiers/dqi/verify_submission_dqi.py PATH/dqi.kmb
python from_host/verifiers/dfthcsos/verify.py PATH/dfthcsos
```

Keep default sizes, shots, and confidence for final measurements. The metric's `verification_shot_divisor` is testing-only and changes confidence and scores. Official split salts remain private; an independent local salt checks robustness without reproducing the leaderboard samples.
