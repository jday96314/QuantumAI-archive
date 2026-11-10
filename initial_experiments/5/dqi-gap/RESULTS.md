# Round 5 gap-aware compact DQI decoder

**Accepted candidate:** `candidates/gap_best_n70_m4095.kmb`, built by `build_best.py`. General quantum dispatch and the complete phase/inverse oracle are implemented. The default native DQI cost is **1435930717449.9397**, a **1.113009x improvement** (10.15349% lower cost). The preferred 1.25x target was **not reached**; arithmetic T is 17.25% above the incumbent, compensated by the reduced repetition factor.

All **three independent 102,400-shot native runs** have zero value errors and zero phase errors, at unchanged delta=2^-128. Results are not pooled for confidence and no favorable salt is selected. The two or three residual events per million remain explicit in the classical and adversarial records; this is an approximate decoder, not a proof of error-free decoding.

| Metric | Round-4 incumbent | Best complete gap-aware candidate |
|---|---:|---:|
| Q | 961 | 1004 |
| Mean executed T | 1,019,865.856484375 | 1,195,802.154687500 |
| Phase / value failures, default 102400 shots | 380 / 0 | 0 / 0 |
| R | 1.6968460732519504 | 1.1912604648115832 |
| DQI cost | 1,598,204,229,384.9265 | 1,435,930,717,449.9397 |

`CANDIDATE.json` is the concise disposition; `results/best_validation_summary.json` contains all independent native results. The remaining sections retain the experimental history, including rejected designs.


Completed implementation and native verification. All outputs are local to this directory. Round-4 and `../dqi-repair` are read-only dependencies. No other active experiment directories are consulted.

Incumbent: Q=961, mean executed Toffolis 1,019,865.856484375, 380 phase failures / 102,400 shots, R=1.6968460732519504, cost=1,598,204,229,384.9265. At that intermediate checkpoint no replacement had been accepted.

The imported `census.cpp`, classical reference, component builders, and checks are copies of the prior repair work. Their outputs are redirected here. The inherited compact Euclid and field arithmetic attribution is retained in `../dqi-repair/RESULTS.md` and the pinned Apache-2.0 source in `../../3/vendor/uncompute`.

First reproduce the exact-host census and local reversible checks, then test shared normal/exception updates before full integration. Every controlled arithmetic gate is charged; exception frequency never discounts static gates.

## Reproduced evidence

`results/host_census_assessment.json` exactly reproduces 744 wrong-decode events and 380 phase failures with the incumbent's original independent masks. The prior classical splice leaves 2. The imported six-moment, k=10 cofactor, and N70/k10 full gap-block native checks reproduce their published Q/T and pass forward, inverse, and compositions.

## New architectures and correctness

Architecture A (`gap_serial.py`, `quantum_control.py`) uses coherent ordinary/wait/quadratic-carry dispatch with shared polynomial tail passes. Waiting retains the two consumed remainder moments; a quadratic quotient finishes in the next slot. A previous cofactor short by one degree is stored multiplied by x and corrected by reversible coefficient rotations. A separate terminal-shift invariant handles immediate completion and the final two-slot boundary. Completion remains common-x padding. Scratch is retained through phase extraction and cleared by the inverse.

Architecture B (`unified_cofactor.py`, `fused_gap_block.py`) replaces the two selected cofactor arithmetic networks with one: beta*P+(alpha+lambda*x)*B and a routed monomial. The short convolution uses k+1 field products instead of 2k. Three local Schur slots share one pair of large tail passes, including incoming/outgoing carry boundary terms. The native component includes gaps in each of the three slots, incoming carries, shifted previous cofactors, completion, and inverse-first compositions.

A third improvement (`onepass_gap.py`) reconstructs consumed remainder words from surviving cofactor/normalizer invariants. It removes the copied head and two extra cofactor passes from a serial slot. Its N22/k3 test passes 256 cases in all four directions/compositions (Q335, T5317 per direction). Full-size integration and independent verification were subsequently completed below.

### Classical census of the carried representation

The new independent C++ implementation is `carried_census.cpp`. It retains the incumbent model and exact EEA/BM cross-check from the original census.

| Stream | Cases | Residual wrong decodes | Residual random-mask phase errors | EEA/BM disagreements |
|---|---:|---:|---:|---:|
| 2026102405 | 1,000,000 | 3 | 0 | 0 |
| 2026102406 | 1,000,000 | 2 | 1 | 0 |
| Original host syndromes | 102,400 | 0 | 0 | 0 |

The first two residual populations include respectively two and one degree-three quotients, plus one adjacent repeated-gap event in each stream. The decoder repairs separated repeated gaps. A zero phase count on three wrong decodes is sampling variability, not exact correctness or an acceptance score. Native full verification remains required.

### Native results so far

- Direct six-moment Mobius shears: all 2/4/6 moment, k1/k3/k10 checks pass, including zero leading moments, inverse and composition. This eliminates unnecessary conjugating polynomial passes.
- Selected serial component: N14/k1, N14/k3, N70/k10, N70/k33, N70/k34 pass all four directions/compositions and clean scratch/phase.
- Unified cofactor component: N14/k3 and N70/k10 pass the same checks.
- Fused component N22/k3: 256 cases pass all four modes, including every local gap position and carries across the block boundary. Initial implementation used the inverse scale in the wrong direction; rejected logs are preserved (`fused_gap_initial.log`). A concurrent check overwrote intermediate circuit names once; later checks use distinct `DQI_GAP_TAG` suffixes.
- Full selected serial N14/M63: unchanged native verifier, 4096 shots, zero value/phase errors. A 1024-shot diagnostic had zero errors but fails the host's fixed confidence threshold; it is not claimed as acceptance.
- Full fused N14/M63: unchanged native verifier, 4096 shots, zero value/phase errors.

### Resource failures and continued optimization

The complete selected serial N70 builder produces Q1039 / static T2,984,206. It is too expensive. Nested controls around already-selected affine padding were removed; the original partial build log is retained. This is an implementation cost failure, not a mathematical obstruction.

The complete fused N70 builder produces Q1057 / static T1,482,764, a substantial reduction but still not competitive. Its 136,661,248 core operations, allocation trace, and complete builder/circuit are preserved. The subsequent integrated candidate uses the proven one-pass slots at later stages and reuses cleared unary completion flags as phase-scan accumulators. Resource-based architecture selection precedes full native sampling; no circuit-hash sample is used to choose a candidate.

## Current reproduction commands

Run here, with `source env.sh`. All generated files remain here.

```bash
./census 102400 2026102405 parameters/host_cases.bin
"$EXP_PYTHON" host_census.py --assess
"$EXP_PYTHON" check_repair_component.py
"$EXP_PYTHON" gap_cofactor.py --k 10
"$EXP_PYTHON" check_gap_block.py --n 70 --k 10
./carried_census 1000000 2026102405
./carried_census 1000000 2026102406
./carried_census 102400 2026102405 parameters/host_cases.bin
DQI_GAP_TAG=_onepass "$EXP_PYTHON" check_gap_serial.py --n 70 --k 22 --onepass
DQI_GAP_TAG=_fused "$EXP_PYTHON" check_fused_gap.py --n 70 --k 31 --count 3
"$EXP_PYTHON" profile_hybrid.py --after 1
"$EXP_PYTHON" profile_hybrid.py --after 100
```

At that intermediate checkpoint no replacement had been accepted. The original incumbent remains unchanged.

## Final implementation and resource selection

The completed builder is `build_best.py`, with `hybrid_decoder.py` integrating the full forward decoder, runtime-mask phase extraction, and exact inverse. Its frozen pre-sampling source/parameter manifest is `parameters/best_candidate_pre_native.json`; the binary/source manifest is `parameters/best_candidate_manifest.json`.

The selected schedule, in degree slots after the opening slot, is:

- Three-slot blocks beginning at k=1,4,7,10.
- Two-slot blocks beginning at k=13,15,17.
- One-pass single slots k=19 through 34.

`profile_schedule.py` exhaustively profiles the implemented one-, two-, and three-slot primitives at every admissible stage and uses dynamic programming to minimize R*Q²*T over their schedules. It enumerates the actual peak-width breakpoints; no arbitrary width cap is imposed. The selection uses the host zero-error confidence floor, not a circuit-hash sample. Two later lifetime improvements reuse an existing history bit as the previous-shift control and temporarily reuse the completion bit as the padding control. They reduce the selected peak from 1006 to **1004**, without changing emitted T or the selected schedule.

Final compiled resources: **Q=1004, static T=1,195,890**, split into 824,740 core and 371,150 scan Toffolis. The binary SHA256 is `8a0675d51edd32a4eed4f5e35e2b3c5e13c83639bdb14142e36c119f7175fde9`. These are measured compilation bounds; the native executed means and accepted scores are recorded below.

### Laurent cofactor and one-pass boundary

`laurent_cofactor.py` removes the conditional current-bank rotations. On a carry, the stored current polynomial has zero constant coefficient, so division by x is exact. A masked alpha word and a Laurent left factor share all ordinary and exceptional products. Evaluation at k nonzero field points uses k field multiplications; the single high coefficient is corrected by a free linear map. The same arithmetic network handles ordinary, waiting, and quadratic-carry charts.

`onepass_gap.py` reconstructs consumed raw remainder words from the output scale and high cofactor coefficients. Its boundary correction is a short polynomial convolution rather than three separate product sweeps. This removes both extra cofactor passes and the copied head from each late slot.

The final Laurent N22/k3 fused component passes 256 cases in forward, inverse, both, and inverse-first modes (Q425, T19,818 per direction). The one-pass N70/k22 component passes the same modes (Q921, T13,051 per direction). Component widths include their public operands and are not full-decoder widths.

### Phase workspace reuse

The scan computes the completion count, XORs the expected unary completion pattern into the stored flags, and uses 24 cleared flags as derivative/previous-value accumulators. Eleven more flags plus one clean qubit hold the explicit leading coefficient during unpadding. The carry word supplies the root accumulator. All rotations and flag XORs are reversed before the decoder inverse.

On a rare non-unary residual path these accumulators may be dirty. They are used only as inputs to evaluation and the phase kernel, which restore their inputs; no dirty bit is measurement-erased as though it were clean. This preserves syndrome and private scratch even on deliberately uncorrected cases.

The full-size shortened-scan diagnostic (`gap_hybrid_scan_test.kmb`, 70 syndromes / 63 mask positions) passes 1024 adversarial cases, every gap stage 0–33, immediate completion, zero syndrome, one-error, weight-35, and random controls, including twice-executed restoration. It is a component/integration diagnostic, not a full benchmark score.

The final full-mask candidate separately passes the same 1024-case adversarial contract, including twice-executed phase/input/scratch restoration. `results/best_adversarial.json` is the independent native report.

### Explicit residual-path preservation

`extract_residual.cpp` regenerates exactly all five wrong-decode events from the two million-case streams, without selecting new samples. Three have a degree-three quotient; two have adjacent degree-two quotients [2,2]. `verify_residual.py` applies 128 independent masks to each:

- All 640 shots restore every syndrome and private scratch.
- Repeating the oracle restores phase and data on every shot.
- Truth-phase errors are 55,63,58,63,57 / 128, respectively, as expected for retained decoding failures.

These failures are preserved explicitly rather than hidden by the ordinary-input tests. The candidate repairs separated repeated gaps. Extending the current routing to degree-three and adjacent-double-gap charts would charge additional arithmetic at every slot for an observed incidence of only 2–3 events per million. No universal resource lower bound for such a correction is claimed.

### Why the 1.25x preference remains unmet

The zero-error host floor at 102400 shots and delta=2^-128 is R=1.1912604648115832. At Q1004, 1.25x improvement requires mean T <= **1,064,751.13**. The current circuit is about 131,000 Toffolis above that budget. This is an implementation limitation, not a mathematical incompatibility.

The remaining cost is visible in `results/profile_parts_k10.json`: the incumbent degree-10 three-slot block costs 34,953 T, while the completed gap-aware block costs 43,470 T. The large tail products cost the same. The added cost is concentrated in cofactor routing (14,049 versus 8,982 T), general local Schur/dispatch work, and incoming/outgoing carry normalization. Blocked updates still copy a seven-word head and use three cofactor passes to expose both old and new polynomial frames. The one-pass reconstruction avoids this late in the decoder, but using it early repeats long tail products. The schedule search prices that tradeoff explicitly.

The actual 1004-qubit peak is recorded in `results/gap_best_n70_m4095_allocations.json`, together with all live windows at the peak and an empty final allocation list. `results/gap_best_n70_m4095_lifetimes.json` records per-stage emitted Toffoli counts and local peaks. The peak occurs in the k=10 Schur head: 840 syndrome bits, 26 history bits, 12 carry bits, two persistent controls, three emitter bits, an 84-bit copied head, one normal-chart control, and three 12-bit field windows. This ledger explains the current width without treating it as a proven minimum.

### Native benchmark validation

`validate_best.py` runs the unchanged host verifier for 102400 shots each with the default hash seed and independent salts `dqi-gap-validation-20261024-a` and `dqi-gap-validation-20261024-b`. The candidate was fixed before all three runs. A partial runtime attempt with two workers was stopped before producing results; all fixed streams were restarted with four workers solely to reduce runtime. No samples, confidence parameters, verifier code, or candidate were changed, and no partial result was used for selection.

All three unchanged native runs passed; their counters are recorded below. Independent validation results remain separate; confidence is not pooled and the most favorable salt is not selected.

### Additional reproduction

```bash
source env.sh
# Rebuild the frozen resource-selected circuit, then run all native checks.
"$EXP_PYTHON" build_best.py
"$EXP_PYTHON" validate_best.py
# Explicit residual cases and reversible-preservation checks.
g++ -O3 -std=c++17 extract_residual.cpp -o extract_residual
./extract_residual
"$EXP_PYTHON" verify_residual.py
```

The initial classical reference log `serial_reference_exceptions.log` used an incorrect exact-EEA stopping condition; `serial_reference_exceptions_v2.log` is the corrected reference. This was a test-oracle error, not evidence of a decoder failure. Earlier rejected component logs and expensive full circuits remain preserved.

## Accepted native benchmark results

| Fixed stream | Shots | Phase / value errors | Mean executed T | R | DQI cost |
|---|---:|---:|---:|---:|---:|
| default | 102400 | 0 / 0 | 1,195,802.154687500 | 1.191260464812 | 1,435,930,717,449.9397 |
| salt_a | 102400 | 0 / 0 | 1,195,802.071328125 | 1.191260464812 | 1,435,930,617,351.2012 |
| salt_b | 102400 | 0 / 0 | 1,195,802.208828125 | 1.191260464812 | 1,435,930,782,462.5222 |

The **default stream** supplies the accepted score and 1.1130092907x improvement. The two salts validate the same frozen binary. All three use effective error bound 0.0008660587301657926 at delta=2^-128.

Final integrity audit: the host verifier source hash, the round-4 incumbent binary hash, and every source in the best-candidate manifest are unchanged. All allocator windows close; the measured allocator peak and native circuit width both equal 1004. No previous candidate or other active experiment directory was modified.
