# Optimization iteration 5

The strongest saved round-five package is **[integration/wave3/submission.zip](../initial_experiments/5/integration/wave3/submission.zip)**, with local default score **9.811059643318059**, compared with round four's **9.367722153134604**. The useful changes were gap-aware DQI decoding, cheaper middle-product arithmetic, a small chemistry normalization improvement, and a chemistry lookup compiler that reduced scored width. No inverse-square-root experiment beat the inherited round-four circuit.

This report summarizes **all 25 Markdown reports** in [brainstorming/results](../initial_experiments/5/brainstorming/results), across three waves, and reconciles their intermediate selections with the actual final integrated archive. The individual reports include overlapping work and historical status text; their gains must not be counted twice. Detailed native reports, certificates, failing controls, and reproduction commands remain in the linked source documents and experiment directories.

## Best saved round-five submission

| Task | Round-four selected cost | Final round-five cost | Selected improvement |
| --- | ---: | ---: | --- |
| Inverse square root | 182,956.5271606425 | 182,956.5271606425 | Retain `zero_head_affine.kmb` |
| DQI | 1,598,204,229,384.9265 | **1,414,354,472,410.0034** | Gap-aware decoder with trimmed middle products |
| Chemistry | 708,547,388.9553173 | **696,942,846.9115301** | Banked physical fit plus structured four-word lookup |

The score is `(2076774.233715/inverse_cost * 24439463583672.93/DQI_cost * 3355599588.315537/chemistry_cost)^(1/3)`. Lower task costs and higher combined scores are better. These are local full-workload scores using supplied baseline costs, not private competition results.

The [final wave-three integration report](../initial_experiments/5/integration/wave3/RESULTS.md), [summary](../initial_experiments/5/integration/wave3/submission.summary.json), [archive audit](../initial_experiments/5/integration/wave3/results/archive_audit.json), and [original-metric audit](../initial_experiments/5/integration/wave3/results/metric_preflight_and_score.json) bind the exact selected bytes to full-shot native evidence. Archive SHA256:

```text
b521c86dfdf830ada6e0f9cbfc88a5e58c762be3a2feabaa50c0479324652fa1
```

The earlier [integration/submission.zip](../initial_experiments/5/integration/submission.zip), score **9.709873458706737**, remains preserved. It combined the original gap decoder with banked chemistry. The chem-audit workspace also packaged a valid **9.4712493866** combination using the trimmed older decoder; that was not the strongest round-five package.

## Wave 1: diagnosis and first architecture probes

### DQI phase extraction

[dqi_phase_extraction_report.md](../initial_experiments/5/brainstorming/results/wave1/dqi_phase_extraction_report.md) separated the locator representation from the cost of applying the phase. It independently established field-basis conversion, the trace-dual mask map, and the distinction between the compact prior-cofactor state and a standard evaluator. Treating that cofactor as the evaluator would be incorrect.

The arbitrary-mask interpolation construction guarantees **70** skipped positions, not the fixed-mask 71st skip sometimes quoted elsewhere; expected active positions are about **4,024.0173**. Classical interpolation has substantial instruction overhead even though it costs no quantum width or Toffolis in this objective.

| Tested phase route | What worked | Why it did not replace the decoder |
| --- | --- | --- |
| Prepared amplitude F | **11-T point body**, runtime masks, root gating, clean scratch, and phase passed native checks. | The compact state did not already contain F; preparation/storage dominated. |
| Separate evaluator and inverse derivative Psi | **23-T bilinear point body** passed. | Half-degree Psi preparation cost **88,930 T** for build/erase and **1,383 live Q**, before missing evaluator conversion. Favorable screen was already **1.708x incumbent cost**. |
| Standard division-based Forney | 80-T body passed. | It requires a different numerator; per-point conversion would add a full inverse lifetime. |
| 87-T cubic point body | Reproduced valid component and older full result. | Full sampled result cost **1.615697e12**, worse than selected **1.598204e12** after repetition penalty. This did not prove a worse underlying phase function. |
| Schoolbook/fast modular F maps | Reversible modular MACs passed, including arbitrary incoming targets. Fast degree-35 MAC **23,352 T**. | Complete in-place F lifetime **93,408 T**, **2,931 Q**; out-of-place still **2,511 Q**. Missing compact adapters only add cost. |
| Fermat polynomial inverse | Explicit compute/copy/uncompute was constructed and accounted for. | Millions of Toffolis and thousands of scratch qubits; uncompetitive. |

The useful lead was direct low-workspace preparation on the existing polynomial bank. The implemented schedules lost; the 11-T body itself remained valuable and later round-six work successfully revisited it.

### DQI exception diagnosis

[dqi_repair_exception_diagnosis_report.md](../initial_experiments/5/brainstorming/results/wave1/dqi_repair_exception_diagnosis_report.md) reproduced the selected decoder's **380 phase errors/102,400 shots** exactly. Independent EEA and Berlekamp–Massey references agreed on every census case. Reconstructing the original host sample found **744 wrong decoded events**, whose original masks explain all 380 phase failures.

The dominant fault was an interior missing leading coefficient requiring a degree-two quotient; synthetic zero-pivot completion incorrectly treated it as stopping. Million-case streams found **7,135/7,483 wrong vectors**, mostly isolated degree-two gaps, with rarer multiple gaps and degree-three quotients.

A reversible degree-(2,1) repair and cofactor realignment passed native arithmetic, inverse, and phase checks. A forced stage-10 full diagnostic fixed **448/448 targeted exception shots** and saved **10,506 T**. It nevertheless failed every zero-/one-error control's value/scratch restoration and ordinary weight controls. It lacked a quantum dispatcher, immediate-completion handling, and final-boundary support, so it was **rejected as a submission**. A classical splice reduced original-mask phase errors from 380 to two; that was diagnostic evidence, not a native replacement score.

### DQI register design

[dqi_hybrid_register_design_report.md](../initial_experiments/5/brainstorming/results/wave1/dqi_hybrid_register_design_report.md) rebuilt the actual selected round-three and round-four tapes byte-for-byte and measured live registers. The compact state already stores quotient/cofactor information inside the original 840 syndrome qubits; it is not simply an extra history bank that can be discarded.

A reversible discrepancy overwrite component passed native phase, scratch, and composition checks, but parking information was necessary. A tested fiber showed why overwriting without retaining information was not reversible. BM history rose toward Q1428, and no complete hybrid below the requested Q1200 design threshold emerged. Work stopped before rewriting a full decoder. This was a limitation of the proposed schedules, not a universal width bound.

### Chemistry quantization-aware fitting

[chem_qfit_quantization_aware_factorization_report.md](../initial_experiments/5/brainstorming/results/wave1/chem_qfit_quantization_aware_factorization_report.md) calibrated against the **actual accepted integer dictionary**, rather than an already-invalid float control. Fixed-word refits reached about **19.5330 mHa**; ridge regularization changed that only slightly. The target's small negative eigenvalues were retained rather than silently clipped.

Methods tried included tensor-feedback freezing of vectors/angles, hardest-first/random orders, exact coupled integer moves within and across vectors, residual-mode whole-row exchanges, legacy block descent, annealing/QUBO controls, leverage pruning, integer column splitting, normalization/offset constraints, and gauge alignment between coefficient refits. Some improved residuals: the accepted-basin discrete control reached **19.5163741 mHa**, and B671 exchange reached **23.7993 mHa**. Smaller two-block dictionaries remained invalid. Matched gauge-alignment ablations showed that alignment matters but did not solve the integer accuracy problem.

| Fully native-tested new fit | Physical result | Cost | Selection |
| --- | --- | ---: | --- |
| Greedy B711/C350 | Passed both accuracy gates | 815,675,394.36 | Rejected on cost |
| Integer B700/C360 | Passed both accuracy gates | 825,280,181.62 | Rejected on cost |
| Constrained B711/C350 | Passed both accuracy gates | 728,405,368.77 | Rejected on cost |

The smaller dictionaries' normalization increases outweighed preparation savings. Full workload shot counts vary with dictionary shape; these were full native workloads, not shortened tests. The round-four **708,547,388.96** selection remained stronger.

### Chemistry multiblock and generator searches

Two linked reports cover distinct finite searches:

- [chem_blocks_multiblock_search_report.md](../initial_experiments/5/brainstorming/results/wave1/chem_blocks_multiblock_search_report.md): exact host resource frontier; R4/B180/C88 random/PCA warm starts, R2/B351/C175 PCA, covariance refits, seven-bit quantization and one eight-bit fallback. Best continuous R2 error **27.153 mHa**; best tested rounded/refitted two-body errors still around **172.525 mHa** or worse. Both accuracy and normalization failed.
- [chem_blocks_jad_generator_resource_search_report.md](../initial_experiments/5/brainstorming/results/wave1/chem_blocks_jad_generator_resource_search_report.md): **28,996** constructively costed configurations, signed difference-bank oracles, physical generator-gauge/JAD controls, R4/B180/C88 and a single R4/B224/C88 eight-bit pivot. Continuous errors **54.578–65.952 mHa** in B180 and **31.855 mHa** in B224; best quantized covariance errors about **273.69/153.94 mHa**. No replacement.

Conditional resource estimates such as 368–489M assumed incumbent normalization and physical feasibility. They were not achieved scores. The signed feature spectrum's necessary rank near 290 screens out some shapes but does not establish feasibility of the remaining shapes. Native float32-vs-pair-space target discrepancies were measured and corrected before final comparisons; failed precision guards are retained in the source evidence.

## Wave 2: gap decoding, arithmetic, and stronger chemistry certificates

### DQI dispatch, core trimming, and gap-aware decoder

| Report | Successful work | Final outcome |
| --- | --- | --- |
| [Claude: exception dispatch](../initial_experiments/5/brainstorming/results/wave2/claude/dqi_dispatch_exception_correction_report.md) | Autonomous classical local-transfer model decoded all 102,400 original cases; million-case streams had one/zero wrong vectors. Matrix interiors, Bennett bridge, and completion-counter components passed. | Matrix/bridge schedule added **115,482 T** even before cleanup and **176,466 T** including another cleanup. Killed on resources; no complete accepted decoder. |
| [Claude: compact core](../initial_experiments/5/brainstorming/results/wave2/claude/dqi_core_compact_euclid_report.md) and [design memo](../initial_experiments/5/brainstorming/results/wave2/claude/DESIGN_MEMO.md) | Trimmed zero rows/columns and retained immutable inverses in eligible contexts: core **649,802 -> 624,206 T**, saving **25,596 T**, same Q961. Twenty-eight native component circuits passed; 120 alternatives and dynamic programming checked schedules. | At this wave, no full candidate or acceptance; the later wave-three report supersedes that stopping status. Block-4 schedules exceeded width, and repeated whole-series running-ratio inversions lost badly. |
| [GPT: gap-aware decoder](../initial_experiments/5/brainstorming/results/wave2/gpt/dqi_gap_aware_decoder_report.md) | Full quantum dispatch, carried representation, Laurent cofactor, boundary handling, and workspace reuse. **Q1004**, static **1,195,890 T**. All three fixed **102,400-shot** runs had zero value/phase errors. | **Accepted cost 1,435,930,717,449.9397**, 10.15349% below round four. Higher T and Q were offset by repetition **1.69685 -> 1.19126**. Preferred 1.25x gain not reached. |

The gap decoder still has explicitly retained rare degree-three/adjacent-double-gap residual cases: independent million-case censuses found three/two wrong vectors. All gates are charged even when exception controls are rarely active. Early scale-direction and reference stopping-condition mistakes were corrected and their rejected logs retained. Native adversarial and twice-executed restoration checks covered the supported dispatch cases and known residuals.

### Chemistry normalization audit

[chem_audit_assumptions_and_margin_report.md](../initial_experiments/5/brainstorming/results/wave2/claude/chem_audit_assumptions_and_margin_report.md) produced a small genuine gain: **708,547,388.9553 -> 708,156,073.7876**, with unchanged T30,504/Q478 and improved effective normalization **48.56733747**. The banked physical fit passed full native verification and became the physical input for later compiler composition.

The report also derived a bound applying to every valid factorization under the fixed verifier's normalization/SOS definitions and physical domain:

```text
lambda_eff > sqrt(3.914898 * 559.1324731485038)
           = 46.78617959252637
```

It is not restricted to one dictionary or gauge. At fixed incumbent T/Q, it implies cost above about **682.185M**, so normalization alone cannot reach 567M. It does not bound other circuits or the combined score.

Further normalization-aware scratch fitting tested R1/B719, R1/B639, R2/B447/B463/B479 and R4/B180 for 20,000 steps, followed by rounded/refitted physical checks. Best float errors were roughly **27.8–45.9 mHa**, and low-Lambda checkpoints remained inaccurate; decoded refits were roughly **128.8–153.9 mHa**. A differentiability/precision check and GPU-memory issue were corrected. The prepared noise/descent/growth branch was **not run** under the experiment's stopping rule.

### Direct discrete chemistry and block-PSD fitting

[chem_discrete_integer_optimization_report.md](../initial_experiments/5/brainstorming/results/wave2/gpt/chem_discrete_integer_optimization_report.md) tried exact integer curvature proposals, simultaneous QUBO moves, subtree annealing, residual-geometry descent, donor/receiver exchanges, rank-aware screens, PSD-geometry updates, and normalization-constrained native one-body matching. Numerical gradient and rank-exchange controls passed. Better words did not become a cheaper feasible model.

Three feasible B711 native fits cost **727.980M, 727.915M and 727.841M**, all worse than 708.547M. A matched B479 direct-search control showed tiny error changes around **22.547 mHa**; annealing accepted real uphill moves but found no passing endpoint. The saved B447 block-linear floor was **23.3047 mHa**. For the original B479 words, unlimited-rank PSD primal/dual residuals were **20.156416/20.156390 mHa**, excluding coefficient-only rescue of that dictionary. Changing its words can beat that floor, so it is not a shape-wide theorem.

The strongest new actual C175 B479 dictionary reached **22.402124 mHa**, with normalization **73.4511**, still invalid and expensive. Relaxed PSD-geometry descent lowered the old floor to **20.067151 mHa** using C479, while actual C175 audits remained around **22.5860 mHa**. Independent row-support/gauge initializers plateaued near **46.3811/235.1558 mHa** for B479/B447. Whitened coefficient L-BFGS reached **22.4178 mHa**; nonconvex rank-175 ADMM was unstable. B639 continuation improved **102.0062 -> 34.0700 -> 30.2548 -> 29.8437 mHa**, also invalid. Whole-atom uphill searches and quenches retained their original best states. These are different failure modes, not evidence that the entire shapes are impossible.

[chem_psd_block_factorization_report.md](../initial_experiments/5/brainstorming/results/wave2/gpt/chem_psd_block_factorization_report.md) separated fixed geometry from covariance rank. Rank continuation, unrestricted PSD/dual solves, projected-trust geometry changes, graph/eigen initializers, integer Newton/trust steps, and residual-directed augmentation were tested.

| Result | Meaning |
| --- | --- |
| R2/B351/C189 continuous error **19.978781 mHa** | Valid continuous construction, conditional 1.4047x resource gain at incumbent normalization; not a native submission |
| Same rounded dictionary union-span floor **140.196114 mHa** | More coefficient columns cannot rescue those decoded vectors |
| Continuous trust endpoints about **22.137/23.865 mHa** for strongest tested prior geometries | Still above accuracy before quantization |
| Exact-word R2/B375 residual pivot full-PSD error **51.459643 mHa** | Best new integer geometry still invalid; explicit dual floor above 51 mHa |
| Actual R2/B375/C192 word screen **53.181669 mHa** | Failed native physical gate despite resource appeal |

All actual ANGLE families for the physical screens were compiled and full resource accounting checked. Invalid physical products were not reported as accepted scores.

### Inverse width and adaptive dyadic designs

| Report | Tested methods | Result |
| --- | --- | --- |
| [Claude: inverse width](../initial_experiments/5/brainstorming/results/wave2/claude/inv_width_inverse_square_root_report.md) | Byte-identical baseline, exact synthesis margin, low-ancilla adders, coefficient delays/curvature release, segmented phasing and measured carry erasure, dirty input borrowing, direct quadratic, bit-serial Horner, and aligned linear tables | Best incumbent schedule **Q88/T1842**, projected cost **183,202.82**, worse than 182,956.53. Other full-design projections roughly **448k–2.02M**; no full acceptance run. Narrow precision variants failed accuracy. |
| [GPT: adaptive dyadic inverse](../initial_experiments/5/brainstorming/results/wave2/gpt/inv_dyadic_inverse_square_root_report.md) | Twenty integer fits, local/shared-coordinate dyadic partitions, mixed linear/quadratic tails, address refinement, coefficient/lifetime sweeps, streamed ANF erasure | Both frozen candidates passed all **200,540,928** shots. Shared Q90 cost **182,956.5389** merely matched the resource point; local Q92 cost **188,762.1693** lost. |

The reusable dyadic component streamed an 11-bit-address ANF bank with **31 rather than 83 scratch wires at the same 83 T**. It improved losing linear/mixed schedules but did not lower the incumbent's peak. Accurate approximation, lower component width, and complete objective improvement were kept distinct.

## Wave 3: accepted compiler/core changes and architecture census

### Accepted trimmed decoder and its later integration

[dqi_core_trimmed_decoder_report.md](../initial_experiments/5/brainstorming/results/wave3/claude/dqi_core_trimmed_decoder_report.md) **supersedes its embedded wave-two stopping text**. A complete Q961/static-T994,358 decoder was built and accepted. A paired-seed comparison retained exactly the same 395 phase failures, width, and repetition factor while reducing mean Toffolis by **25,596**. Default cost was **1.5599009371e12**, with 382 phase failures; independent cost **1.5571895502e12**, with 379.

This older trimmed architecture still loses to the gap-aware decoder. The final round-five integration instead **ported its trimmed middle-product emitter into the gap decoder**, retaining the gap schedule/dispatch and original inverse. It saved **23,352 static Toffolis** at Q1004: complete count **1,172,538** (801,388 core + 371,150 scan).

The final gap-trimmed default run had **zero value errors and one phase error/102,400 shots**, cost **1,414,354,472,410.0034**. Paired tests and a full original-decoder replay on the new default stream reproduced the same error counts and repetition factor, with exactly 23,352 more Toffolis in the old circuit. Thus the gain is arithmetic, not favorable sampling. See [acceptance](../initial_experiments/5/integration/wave3/results/gap_trimmed_acceptance.json) and [controlled default-stream comparison](../initial_experiments/5/integration/wave3/results/new_default_stream_comparison.json).

### Dispatch and decoder census routes that did not replace it

| Report | What was tried and learned | Final disposition |
| --- | --- | --- |
| [Claude: fixed-schedule dispatch](../initial_experiments/5/brainstorming/results/wave3/claude/dqi_dispatch_exception_correction_report.md) | Classical model passed exact host, two million-case streams and all-weight controls with zero wrong vectors. Shared cofactors, carry-finish slots, routed Schur heads, conditional moments and merged bank routing passed native components. | Q973 is an embedding/plan, not full decoder width. Partial extra **326,622 T**, already over the revised +300k limit, omitted more work. Separate carry-chart moments were the main blocker. |
| [Claude: decoder census](../initial_experiments/5/brainstorming/results/wave3/claude/dqi_census_decoder_architecture_report.md) | Compared eight architectures: direct amplitude F', lazy BM handoff, resident sigma/omega, Hankel elimination, wide Psi, root lists, transposed residues and Wronskians. Ten native circuits, **27,136** final-suite shots passed. | Resident update required Q1707; direct Frobenius F' route required Q11751 and millions of T. A low-memory quadratic reciprocal-product overwrite worked at **960 T, 69 extra Q**, but its degree-35 extension remained unresolved. No 2.5x gain. |
| [GPT: phase-aware decoder](../initial_experiments/5/brainstorming/results/wave3/gpt/dqi_core_phase_aware_decoder_report.md) | Derived Riccati tangent/evaluator and general-quotient Wronskian invariants. Built in-place C-fraction small full oracles, degree-35 H preparation, and **46-T/66-Q** inverse-only point body with **61,440** kernel shots. | C-fraction and padded-H independent 20k censuses had **253/291** and **81/78** wrong vectors. Conditional zero-error gains only **1.06361x/1.03680x**; chart errors erased them. Serial functional >=11,833,500 T lost. No full-size acceptance. |

The phase-aware report corrected initial censuses that used amplitudes instead of the host's squared-amplitude weight law. Those rejected censuses do not support its final conclusions. Degree-35 H build/erase was **100,080 T at Q897**; useful component evidence did not establish a complete supported decoder. Blocked extraction, direct (B,H) updates, common moment maps, and degree-35 low-memory inversion remained research leads.

### Chemistry pruning audit and complete intermediate archive

[chem_audit_pruning_and_submission_report.md](../initial_experiments/5/brainstorming/results/wave3/claude/chem_audit_pruning_and_submission_report.md) continues the banked-normalization report. It tested exact truth-table subtree pruning and address ordering: every one of eleven oracles saved **zero T**, so no compiler replacement followed from that probe.

Its complete archive used banked chemistry and the accepted older trimmed decoder and scored **9.4712493866** in a full host-metric run. DQI's cost there, about **1.54722e12**, reflects that run's sample; it does not supersede the stronger gap decoder. The native inverse cost also differs slightly through measurement statistics despite unchanged bytes. The newer integration uses component evidence and controlled comparisons explicitly.

### Content-aware chemistry QROM and physical/compiler composition

[chem_qrom_content_aware_compilation_report.md](../initial_experiments/5/brainstorming/results/wave3/gpt/chem_qrom_content_aware_compilation_report.md) retained R1/B719/C350 physical accuracy and tested care-prefix erasure, parity pairing, direct streamed ANF, four-row nonlinear banks, local linear/quadratic decoders, batching, row/address permutations, and exact angle gauges. All-address tests included arbitrary incoming output words; they caught and rejected a compressed-bank variant that reintroduced a redundant high payload bit.

| Compiler family | Total task cost | Outcome |
| --- | ---: | --- |
| Baseline J7 | 708,547,388.96 | Reference |
| Care-prefix two-row | 704,970,274.55 | Full native valid gain |
| Parity-paired two-row | 704,737,994.39 | Full native valid gain |
| Smaller J two-row / direct ANF | About 732M–863M | Correct component families, worse objective |
| Local linear four-row | 702,269,799.06 | Full native valid gain |
| Quadratic four-row | 698,078,065.52 | Full native valid gain |
| Frozen hybrid four-row | **697,327,965.84** | Two full native passes; selected compiler |

The winner uses **19 J4 oracles**, T31,608 rather than T30,504, but Q454 rather than Q478. The lower persistent output width wins despite the gate increase. Each oracle fits the 128-Q cap. Rotations, state preparation, and shared temporary workspace were included.

[chem_census_factorization_and_compilation_report.md](../initial_experiments/5/brainstorming/results/wave3/claude/chem_census_factorization_and_compilation_report.md) composed that compiler with the **already banked lower-normalization fit**, reaching **696,942,846.9115301**, with two **50,880-shot** native passes. This is the final round-five chemistry. The 697.328M and 696.943M gains overlap; they are not independent reductions to multiply again.

That census also ranked and tested integer pruning, residual-directed growth/reinsertion, rank reduction, one-body/precision simplification, host threshold changes, and physical gauge continuation. B703 passed accuracy but cost **807.507M** due to normalization **56.3829**. Growth B511 failed at **240.078 mHa**; growth B671 failed at **34.985 mHa**. Relaxed span/PSD bounds excluded weight-only rescue for the frozen growth dictionaries. C340 could reach **19.910161 mHa**, but C320 and below were excluded for that fixed dictionary; arbitrary future geometries were not excluded. Uniform one-body precision reduction failed at **123.108 mHa**. Normalization-only routes could not supply the target 1.25x reduction.

### Shared-tree topology and grid experiments

[chem_tree_topology_and_quantization_report.md](../initial_experiments/5/brainstorming/results/wave3/gpt/chem_tree_topology_and_quantization_report.md) established legal ordered rooted-tree semantics and screened **20 trees / 17 distinct amplitude lattices** on four controls. It tried orbital/sibling/reparent topology search, random balanced controls, local tree changes, legal half-step codebooks, exact integer repairs and coefficient-only continuation.

Topology improved B351 repair to **54.288581 mHa** versus **58.690849** with the incumbent tree, but both failed. New-tree B479 reached **22.500832 mHa**; its unlimited-rank PSD dual floor **20.087505 mHa** excludes coefficient-only rescue of those words. The shifted-grid B351 candidate failed accuracy and widened resources to Q778 despite using brot9/J6. Full native checks rejected the selected new geometries. These are fixed-word/topology outcomes, not an impossibility result for all trees or multiblock models.

### Inverse architecture census and phase-oriented candidates

[inv_census_inverse_architecture_report.md](../initial_experiments/5/brainstorming/results/wave3/claude/inv_census_inverse_architecture_report.md) compared eight architecture families, then built input-borrowed headroom, circular-result wrap correction, and exponent-Walsh selection. The first two projected **195,644.62/190,469.06**, worse than 182,956.53. Walsh reduced width/gates to Q88/T1639 but required **307 rotations**, raising projected cost to **354,397.57**.

Thirty exhaustive integer pipelines found a useful narrow curvature-square model: keep the linear coordinate exact and round only inside the quadratic term. A q5 variant passed with worst error/tolerance **0.8046964**, and a q4 alternative passed more accurately; q6 failed. No native complete fused-square circuit was built in this round. An injective exponent packing model also passed all-input checks but lacked a costed reversible normalizer; no physical five-qubit saving was claimed.

[inv_phase_inverse_square_root_report.md](../initial_experiments/5/brainstorming/results/wave3/gpt/inv_phase_inverse_square_root_report.md) built **two complete, byte-reproducible candidates**, each passing **200,540,928** native shots:

| Complete inverse candidate | Q | Static T | Rotations | Full task cost |
| --- | ---: | ---: | ---: | ---: |
| Inherited selected quadratic | 90 | 1,793 | 34 | **182,956.53** |
| Biased small product bank | 84 | 2,293 | 63 | 231,373.44 |
| Direct signed-row phasing | 80 | 2,434 | 201 | 318,150.29 |

Direct phasing avoided a completed result register but paid more for lookup, exponent routing, coefficient replay, measurement-phase repair, and rotation synthesis. Coarser independent row rounding failed tolerance even at two-unit spacing. Narrower width alone did not imply a cheaper phase circuit. Both candidates remained useful valid controls and were rejected on cost.

## Evidence interpretation and coverage

Every source report in the requested results tree is linked above: six wave-one reports, eight wave-two experiment reports plus the core design memo, and ten wave-three reports. Their underlying work is summarized by method and outcome, including stopping decisions and proposed branches that were never run. Saved native PASS results, native failures, classical censuses, numerical certificates, constructive resource profiles, and unbuilt targets are distinguished throughout.

The final integration's [ten-report decision review](../initial_experiments/5/integration/wave3/results/report_selection_review.json) records the wave-three choices and source hashes. It confirms that the old trimmed decoder itself was not selected, the useful middle-product emitter was ported to the gap decoder, the two chemistry compiler reports overlap, and neither inverse report supplied a cheaper accepted circuit.

For the final round-five package, DQI had full default/independent/paired runs plus adversarial and known-residual restoration checks; chemistry had fresh full verification and matching source passes; the unchanged inverse reused its full **200,540,928-shot** evidence. ZIP member hashes and the canonical chemistry digest bind reused reports to the delivered bytes. No test workload, confidence rule, tolerance, or task-cost formula was changed.

Fixed-dictionary bounds explain why more fitting of the same words often could not work. They do not rule out other integer geometries. Continuous accuracy alone repeatedly failed after Givens-word decoding. Partial gate counts also missed decisive storage, reversal, conversion, or synthesis costs. The successful changes improved a complete objective with supported evidence, even when the individual gain was smaller than the original research target.

The later round-six best is documented separately in [iteration 6](../initial_experiments/6/REPORT.md). The documentation copy of this report is [optimization_iteration_5.md](optimization_iteration_5.md).
