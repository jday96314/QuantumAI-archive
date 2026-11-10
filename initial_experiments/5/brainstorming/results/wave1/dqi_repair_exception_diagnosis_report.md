# Round 5 compact DQI exception repair

**Final status:** diagnosed and locally repaired the dominant degree-gap exception; no globally dispatched replacement was accepted. The full forced-branch circuit passes targeted exceptions and fails ordinary-input controls. Preserve the round-4 incumbent.

Scope: round-4 Codex selected `cubic_trace_product_n70_m4095.kmb` only. All prior assets and host verifier are read-only. No other round-5 experiment directories were consulted.

Baseline (retained, not replaced): SHA256 `4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0`; Q=961; executed mean Toffolis=1019865.856484375; R=1.6968460732519504; 0 value / 380 phase failures in 102400 shots; cost=1598204229384.9265. Selected builder: `../../4/codex/build_cubic_compact_decoder.py`, core `compact_product_decoder.py`, parameters N70/block3/serial_after22/linear_padding/virtual_zero/virtual_lead/low_inverse/lazy_flags. Selected full evidence: `../../4/codex/results/cubic_trace_product_n70_m4095_verify.json`.

Read round-4 consolidated and Codex reports and host `verify_submission_dqi.py`. Contract uses syndromes starting at exponent 1, native GF(4096), independent uniform 12-bit classical masks, eigenvector-squared weight distribution. Value failures include syndrome restoration and private scratch. Phase incidence must be distinguished from decoded-event incidence. Host confidence is unchanged, delta=2^-128.

Baseline fresh full reproduction started before any builder editing:
```
source env.sh
"$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py ../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb --max-workers 4 > results/baseline_reproduction.json
```

Hypothesis under investigation: synthetic zero-pivot completion incorrectly treats an interior degree gap as stopping. Prior classical diagnostics skip such cases; that is evidence of a coverage gap, not a verified attribution of the selected circuit failures.

Attribution: the round-4 decoder adapts Apache-2.0 andigu/uncompute, pinned revision `64b5a13aa249f1360846657392e6c66817b26d0f`; original code and license are retained at `../../3/vendor/uncompute`. The compact Euclid architecture is round-4 Codex code. External repository URL was attempted with the web tool, which could not fetch it; the pinned local public source is available read-only.

## Verified diagnosis

Fresh full baseline reproduction is byte/hash consistent and exactly reproduces **all** published counters, including 380 phase failures. `results/baseline_reproduction.json` is the new native output.

`census.cpp` is an independent C++ implementation in the host's native field (modulus 0x1053). Exact polynomial EEA is checked against a separately implemented Berlekamp–Massey decoder on **every case**. The compact circuit model separately implements blocked ratio/quotient conversion, its affine synthetic-pivot correction, late serial arming, unary flag XOR counting, fixed-bank rotations and virtual leading coefficients. It is not imported from the quantum builder. Tables and eigenvector-squared probabilities come from the read-only host helpers. Random positions use rejection sampling without replacement; values and masks are uniform. Independent RNG streams keep masks independent of errors. These million-case streams have the host's distribution, not NumPy's identical sample sequence.

| Census | Cases | Wrong-decode events | Phase failures under independent random masks | EEA/BM disagreements |
| --- | ---: | ---: | ---: | ---: |
| Seed 2026102405 | 1,000,000 | 7,135 | 3,595 | 0 |
| Seed 2026102406 | 1,000,000 | 7,483 | 3,717 | 0 |

All exceptional paths have wrong locator roots, rather than an isolated Forney parity error. The first divergence is a missing expected leading remainder coefficient **before valid completion**. At this point the exact EEA requires a degree-two quotient. The compact algorithm instead arms the zero pivot to one and applies its degree-one synthetic completion correction. Retained lower remainder/ratio coefficients are still nonzero. Subsequent state is therefore no longer the correct RS cofactor state. Gaps occur throughout stages 0–33, with no single stopping boundary accounting for the errors. The field inversion convention and phase kernel are not the first divergence.

In the first million cases, 7,111 events are a single isolated degree-two quotient. Twenty-two cases contain multiple gaps; two contain a degree-three quotient (two consecutive missing leading coefficients). In the independent million, corresponding counts are 7,450 / 32 / 1. There is a dominant reproducible exception, but it requires arbitrary-stage dispatch, not just repairing the first pivot.

`host_census.py` additionally reconstructs **exactly** all 102,400 original verifier shots, using the circuit-hash-derived base seed, SeedSequence batching, full independent mask arrays and original RNG ordering. Every generated syndrome, mask and phase is compared with unmodified `sample_shots`. The independent census finds **744 decoded-error events** in these original samples. Applying the original masks to the model's decoded values gives **exactly 380 phase failures**, explaining the original native count. The `phase_errors` field in `census_host_exact_final.json` uses an independent C++ mask stream; `host_census_assessment.json` / `host_census_assessment_final.log` supply the actual original-mask comparison. This distinction is deliberate.

Native model controls:

- The first model was rejected: it disagreed on 29/256 phases because late serial synthetic pivots were incorrectly modeled as blocked ones. `native_crosscheck_rejected_v1.log` preserves this failure.
- After correcting that distinction, the 70-exception/256-mask native check has zero disagreement (`native_crosscheck_v2.json`).
- `native_adversarial.json`: **1,024 native cases**, including 768 distinct exceptional inputs and 64 each zero-syndrome, one-error, weight-35 and random-weight controls, all with fresh independent masks. Zero model phase disagreements; complete syndrome restoration and clean scratch. Repeating the complete oracle restores phase as well as inputs and scratch.
- Full compact field-state comparisons against exact EEA also pass 1,024 generic cases per million-case census.

## Correction implemented and its current scope

The new `gap_cofactor.py` and `gap_block.py` implement a **degree-(2,1) realignment**. A zero local moment followed by a nonzero moment must consume a degree-two quotient, then a degree-one quotient. After these three degree slots, the ordinary consecutive-degree compact layout is restored. One exceptional coordinate is retained in the six consumed remainder words, rather than adding a separate field-word history bank. An independent classical inverse and full polynomial EEA comparison are in `classical_repair.py`.

For local moments `(0,d,a,b,c,f)`, write

```
u1 = a/d
u0 = b/d + u1^2
beta = d
gamma = c/d + u1^3
alpha = (f/d + u0^2 + u0*u1^2 + u1^4)/gamma
B2 = (x^2 + u1*x + u0)*B + beta*Bprev
B3 = (x + alpha)*B2 + gamma*B
pinv_new = pinv/(beta*gamma)
```

The new block reconstructs the retained high remainder by the exact polynomial product/division identity, uses three cofactor passes like the selected block, and aligns the resulting physical coefficient banks by reversible word permutations. Every arithmetic operation is emitted and charged regardless of input rarity. There is no assumption that a quantum-controlled fallback is free on ordinary shots.

The branch primitive requires `d != 0` and `gamma != 0`. Zero `gamma` can mean genuine completion immediately after the degree-two quotient. The classical splice handles that case separately, and handles the final two-slot boundary. Those two cases are **not** implemented by the current quantum branch. The branch also lacks a global quantum dispatcher that preserves the normal path and handles a gap at any stage. Consequently these components are research artifacts, not accepted replacement circuits.

| New component | Native checks | Total component Q | Static Toffolis per direction |
| --- | --- | ---: | ---: |
| Six-moment reversible Padé reference | 4,096 each forward/inverse/composition | 126 | 2,094 |
| New cofactor realignment, k=1 | 4,096 each direction and both composition orders | 123 | 540 |
| New cofactor realignment, k=10 | same | 339 | 2,430 |
| New cofactor realignment, k=32 | same | 867 | 7,050 |
| New full compact gap block, N22/k3 | 512 each direction and both compositions | 363 | 9,750 |
| New full compact gap block, N70/k10 | same | 939 | 29,700 |

These are complete component widths, including their public operands, and **not** a complete decoder width. Counts are gate counts of compiled circuits, not a DQI score. Public coefficients, the retained exceptional coordinate, zero scratch and exact phase are checked. Half the cofactor tests use arbitrary incoming exceptional-coordinate values. Inverse-first tests use the valid inverse domain; arbitrary singular pivot states are not claimed to be supported. Initial cofactor implementations had top-coefficient ordering/index errors; their failing logs and JSON are retained. The passing version is covered by the final `gap_cofactor_k*.json` results.

The six-moment reference is reused round-4 Codex code, with builds redirected locally; the degree-(2,1) cofactor update and complete gap block are new here. The local gap block is independently tested using native-field polynomial formulas, with field-basis conversion only at the quantum boundary.

A classical splice of the new branch into the compact model leaves **24 wrong decodes in the first million and 33 in the independent million**, exactly the multiple-gap and degree-three events above. This includes the separate immediate-completion and final-two-slot logic. These are **classical diagnostics**, not native full-decoder acceptance or measured R. No replacement score is inferred from them.

## Full diagnostic circuit and acceptance boundary

`build_gap10_diagnostic.py` integrates the new branch into the complete selected builder at a fixed stage 10, including the original prefix, scan kernel and full inverse. It uses reversible word permutations to match the selected builder's physical coefficient layout. This forced branch is deliberately domain-limited: it does not dispatch based on a quantum exception flag. `diagnose_gap10.py` compares it with the selected decoder using the same explicit diagnostic seed 2026102405 and the same syndrome/mask cases. It includes ordinary zero, one-error and weight-35 controls to prevent mistaking success on targeted exceptions for a valid replacement.

No replacement is accepted unless it supports the entire input contract and passes the full native distribution under independent salts. The original selected circuit and combined reference 9.3677221531 remain preserved. Chemistry and inverse assets were not modified.

## Native controls and full-cost accounting

All three baseline controls use all 102,400 shots, unchanged confidence and official scoring. They are the **same circuit**, not competing candidates; the lower independent-salt number is not selected as an improvement.

| Baseline control | Q | Executed mean Toffolis | Phase / value failures | R | DQI cost |
| --- | ---: | ---: | ---: | ---: | ---: |
| Original seed reproduction | 961 | 1,019,865.856484375 | 380 / 0 | 1.6968460732519504 | 1,598,204,229,384.9265 |
| Explicit seed 2026102405 | 961 | 1,019,866.334296875 | 395 / 0 | 1.7115755735216907 | 1,612,078,222,772.0327 |
| Independent salt dqi-repair-control-20261024-a | 961 | 1,019,866.09796875 | 359 / 0 | 1.6760899815248063 | 1,578,655,114,214.934 |

`resource_budget.json` applies the actual host Clopper–Pearson function. At unchanged Q/T, zero observed errors would give R=1.1912604648 and cost=1.122009558305e12 (1.4244x), while eight errors would give R=1.2230139976 and cost=1.151917179973e12 (1.3874x). These are **hypothetical budgets**, not measured replacement results. At Q980 and eight failures, 1.3x improvement permits only about 1,046,659 mean Toffolis. A 1.5x improvement needs additional reduction of the Q²T product even after eliminating all errors.

The concrete next architecture is a reversible dispatcher for the retained-zero degree-(2,1) block, with distinct handling for immediate completion and the final two-slot boundary. It must detect the first gap without overwriting the information needed to clear its flag, route either polynomial bank layout into the selected downstream frames, and charge all controlled arithmetic on every shot. Simply replacing each zero by one, disabling the guard, repairing only stage zero, or charging a complete fallback by exception frequency does not solve this problem.

## Reproduction

Run from this directory. The Python environment and all prior modules are read-only dependencies. `translate_prefix` is a local copy of the pinned round-4 executable because this shared filesystem did not support a symlink.

```bash
source env.sh
# Native baseline and independent control (full 102400 shots by default).
"$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py ../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb --max-workers 4
"$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py ../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb --max-workers 4 --seed 2026102405
"$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py ../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb --max-workers 4 --salt dqi-repair-control-20261024-a
# The checked host field tables/probabilities are saved in parameters/.
g++ -O3 -std=c++17 census.cpp -o census
./census 1000000 2026102405
./census 1000000 2026102406
# Exact original-shot reconstruction; census writes exception_cases.jsonl.
"$EXP_PYTHON" host_census.py
./census 102400 2026102405 parameters/host_cases.bin
"$EXP_PYTHON" host_census.py --assess
# Native simulator controls and repair components.
"$EXP_PYTHON" native_adversarial.py
"$EXP_PYTHON" classical_repair.py
"$EXP_PYTHON" check_repair_component.py
"$EXP_PYTHON" gap_cofactor.py --k 1
"$EXP_PYTHON" gap_cofactor.py --k 10
"$EXP_PYTHON" gap_cofactor.py --k 32
"$EXP_PYTHON" check_gap_block.py --n 22 --k 3 --shots 512
"$EXP_PYTHON" check_gap_block.py --n 70 --k 10 --shots 512
# Domain-limited full diagnostic, not a submission or cost win.
"$EXP_PYTHON" build_gap10_diagnostic.py
"$EXP_PYTHON" diagnose_gap10.py
```

## Final disposition

**No replacement accepted; the selected round-4 circuit remains the incumbent.** The exception diagnosis is complete, and the new reversible branch is validated through the complete original syndrome/phase contract on targeted inputs. General quantum dispatch and its cost remain unresolved. This experiment does not meet the strong or stretch success criterion.

The full forced-stage diagnostic is saved as `candidates/gap10_diagnostic.kmb`, SHA256 `279adee24142ae6825ff39038dd48602017985a79f91de22e62b6566b3151a05`. Its builder and source manifest are retained. On the same explicit seed 2026102405 and exactly the same 512 syndrome/mask cases:

| Full circuit | Q | Executed mean Toffolis on diagnostic set | Stage-10 exception phase failures | Stage-10 exception value failures |
| --- | ---: | ---: | ---: | ---: |
| Selected baseline | 961 | 1,019,862.90625 | 217 / 448 | 0 / 448 |
| Forced stage-10 correction | 961 | 1,009,356.90625 | **0 / 448** | **0 / 448** |

The 448 shots use 20 distinct genuine sampled exception inputs with independent random masks. Both circuits also pass twice-executed input/scratch/phase restoration on those targeted inputs. The corrected full circuit has static T=1,009,448 versus baseline T=1,019,954. Thus the new degree-(2,1) arithmetic branch itself reduces 10,506 gates across forward/inverse at this stage; **this does not price its required quantum dispatch**.

The forced correction is **rejected as a standalone replacement** by the independent controls: all 16 zero-syndrome and all 16 one-error shots have value/scratch failures; weight-12 controls have 9/16 phase failures and weight-35 controls have 8/16 phase failures. The baseline has zero failures on these controls. Its R and end-to-end cost are therefore **not accepted / not reported**: a failure rate measured on the restricted exception set cannot be substituted for native sampling. In accordance with the requested order of checks, this correctness failure prevents proceeding to full native acceptance scoring for this diagnostic circuit. `results/gap10_native_diagnostic.json` contains all groups and both circuits' measured counts.

The classical splice, including its additional completion and final-boundary cases, changes the original masks' **380 phase failures to 2** (`host_census_assessment.json`). Across the million-case streams it produces 14 and 18 random-mask phase failures, respectively. This is strong classical evidence for the algorithmic repair, but is explicitly **not a verified quantum result**. Immediate completion occurs inside the proposed repair in 232 and 255 cases per million, respectively; the quantum dispatch must support it rather than treating every gap as a nonzero degree-(2,1) block. See the final census JSONs for exact counts.

The concrete alternative to the rejected forced branch is an in-place quantum-selected transfer between the regular three-degree block and this exceptional degree-(2,1) block, retaining the zero coordinate in the consumed head. It must also select the degree-two-then-stop case, handle the last two degree slots, and reverse those selections before clearing history. The remaining work is this dispatcher and its native resource/acceptance checks, not another field multiplier or mask-kernel optimization. No lucky salt or restricted-shot result was used to replace the incumbent.
