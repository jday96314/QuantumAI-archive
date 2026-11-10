# Optimization iteration 6

As of October 25, 2026, the best combined local submission is **15.222738605389504**, saved at [integration/submission.zip](../initial_experiments/6/integration/submission.zip). This improves agent14's **14.064304039685362** by **8.2367%**. Agent20 contributes its independently verified public inverse-square-root component and its paired chemistry compiler; agent14 contributes its substantially cheaper DQI decoder. The existing independent-sample reports for these same components give **15.22273646664756**.

The integration is complete. Its [summary](../initial_experiments/6/integration/submission.summary.json), [archive audit](../initial_experiments/6/integration/results/archive_audit.json), [original-metric audit](../initial_experiments/6/integration/results/metric_preflight_and_score.json), [independent-evidence audit](../initial_experiments/6/integration/results/independent_audit.json), and [member/source manifest](../initial_experiments/6/integration/results/package_manifest.json) all record PASS. These are local benchmark scores, not a new competition upload or a private-split evaluation.

## Component comparison and integration decision

Lower task costs are better. The combined score is the geometric mean of the three baseline-cost/task-cost ratios, so task costs must not be added.

| Component | Agent14 selected cost | Agent20 selected cost | Integrated choice | Reduction from agent14 |
| --- | ---: | ---: | --- | ---: |
| Inverse square root | 170,242.24959913077 | **138,095.07008539976** | Agent20: Phuc Thai EXP225a | **18.8832%** |
| DQI | **515,978,772,816.4064** | 1,041,180,174,625.8771 | Agent14: raw continued-fraction destructive quotient | 0% |
| Chemistry | 696,942,846.9115301 | **677,583,323.3862098** | Agent20: paired structured ANGLE lookup | **2.7778%** |

Agent20 also has a fully verified recovered-gap decoder costing **1,023,349,162,559.9081**, awaiting packaging within that workspace. It improves agent20's selected decoder but still costs approximately twice agent14's decoder. It is therefore excluded from this integration. There is no need to combine incompatible decoder internals to obtain the demonstrated improvement: the three benchmark components can be selected independently.

| Submission | Local default score | Disposition |
| --- | ---: | --- |
| Agent14 [submission.zip](../initial_experiments/6/agent14/submission.zip) | 14.064304039685362 | Previous strongest complete submission in these workspaces |
| Agent20 [best/submission.zip](../initial_experiments/6/agent20/best/submission.zip) | 12.046515352810363 | Its actual packaged best, despite older log claims |
| Agent20 with recovered-gap decoder | 12.11607972269253 | Component combination verified; not packaged there |
| **[Integrated submission](../initial_experiments/6/integration/submission.zip)** | **15.222738605389504** | **Best packaged combination reviewed here** |

## What was checked for the integrated submission

The integration selects existing circuit bytes and makes no arithmetic or approximation changes. [package_submission.py](../initial_experiments/6/integration/package_submission.py) checks every default and independent report for PASS, the prescribed full shot count, positive finite cost, and the selected component's digest. It also checks that the chemistry NPZ is byte-identical to agent14's selected physical factorization. Only the chemistry ANGLE circuits change.

The ZIP was extracted and reopened; every member and extracted file was hashed. The canonical chemistry digest was recomputed from the extracted circuits, physical arrays, and tree. The unmodified host metric successfully read the inverse and DQI circuits, all **19 chemistry circuits**, and the chemistry NPZ, and reproduced the score. The metric and verifier files were left unchanged. Attribution and the retained Apache-2.0 license are included in the archive; the [NOTICE](../initial_experiments/6/integration/NOTICE.md) identifies the actual selected DQI implementation rather than agent20's superseded public decoder.

| Selected component | Qubits | Static maximum Toffolis / total chemistry Toffolis | Full shots in each saved run | Native result |
| --- | ---: | ---: | ---: | --- |
| EXP225a inverse | 82 | 1,377; 44 rotations | 200,540,928 | Default and independent PASS |
| Agent14 DQI | 892 | 534,575 | 102,400 | Each run: five value errors, zero phase errors; PASS with repetition penalty |
| Paired chemistry | 454 scored; each oracle within 128 | 30,730 | 50,880 | Default and independent PASS |

Existing full native runs were **reused after digest binding**, not rerun during integration. The independent results remain separate; samples were not pooled or selected for a favorable score. DQI is an approximate decoder with rare unsupported extension poles. Its default repetition factor is **1.2131052140444296**. The inverse's worst error/tolerance ratio is **0.9996626563**; chemistry errors are **11.995373/19.989973 mHa**, below the independent 20 mHa limits. These margins belong to the selected verified components and are not claims of exact all-input decoding.

Archive: **56,655,634 bytes**. SHA256:

```text
c1b38fc6782a81f5991ee17de8e9c6aa4a2629d29348603f14bbbd82fe74f834
```

To reproduce packaging and audit from the repository root:

```bash
source initial_experiments/6/agent20/chem_env.sh
"$EXP_PYTHON" initial_experiments/6/integration/package_submission.py
"$EXP_CHEM_PYTHON" initial_experiments/6/agent20/audit_best_metric.py initial_experiments/6/integration
```

The script writes only the integration package and evidence. It does not replace either agent's selection or milestones.

## Agent14: DQI experiments

The principal success was moving from compact Euclidean/gap-aware decoding to quadratic-extension Schur extraction, then simplifying its representation into raw continued fractions. Interpolation uses irreducible quadratics with nodes outside GF(4096). The final raw recurrence is `R_i = a_i + b_i*x + q_i/R_(i+1)`; rational synthesis updates `(N,B) <- ((a_i+b_i*x)N+q_i*B,N)`. This removes much of the normalization, affine-frame, and node-swap machinery of earlier versions.

The [final optimization report](../initial_experiments/6/agent14/OPTIMIZATION_REPORT.md) and [chronological results](../initial_experiments/6/agent14/RESULTS.md) record the following experiment families, including intermediate designs later superseded.

| Experiment | What worked | What failed, lost, or remained incomplete |
| --- | --- | --- |
| Dynamic programming for triangular/Toeplitz maps | Zero-row/column trimming plus searching recursive splits reduced actual middle-product counts. Sixteen native direction/composition circuits passed. An early complete improvement was packaged. | This retained the inherited gap-aware architecture and was later superseded. |
| GF(4096) cubing through GF(16) evaluation/interpolation | Exact **28-AND** cube versus generic 42-AND multiplication; exhaustive field checks, arbitrary-target and controlled native checks passed. Source masking made the controlled version useful in an accepted early decoder. | Naively controlling the linear corrections increased the core to 800,940 T; rejected. The specialized primitive does not replace general multiplication. |
| General multiplier and joint-cover searches | GF8 quartic seed independently reproduced a different valid 42-AND multiplier. Tensor and cover checks passed. | Cap-2/4/6 walks retained 42 ANDs; a common-operand GF8 pair retained 12 and tested GF16 joint cover retained 16. No general multiplier saving or optimality proof. |
| Blocked ordinary continued fractions | Generic polynomial identities passed, and larger blocks lowered tail-only counts from 594,830 T at block 2 to 188,055 T at block 10. | Counts omitted extraction, transfer matrices, determinant work, synthesis, reversal, scan, and zero charts. No competitive complete decoder was established. |
| Initial quadratic-extension Schur extraction | Exact extension inversion, pivot handling, extraction, and independent natural-law censuses validated the representation. Initial extension inverse: 300 T; full extraction: 194,005 T per direction, Q925. | Extraction alone was not a submission; a 100,000-case model retained five unsupported extension poles. |
| Complete Schur decoder without pivot repair | A complete native phase circuit was built at Q949/T885,278. | **756 value errors and one phase error** made its cost **1.6412e12**, worse than the starting gap decoder. |
| Coherent adjacent-node repair | Quantum node descriptors, controlled constant scaling, structured dynamic folds, ghost-phase repayment, synthesis, and full oracle checks passed. | Storing node flags, controlled frame operations, and the 80-T point scan were expensive. A dropped test-header gate was caught and corrected before acceptance. |
| Fixed-node scalar/gap chart | Classical algebra and native extraction handled scalar pivots without separate node swaps; full extraction saved 86,634 T against repaired extraction. A two-node mixed-degree synthesis component passed through tail degree 33. | The standalone mixed-degree component did not supply a complete dispatcher. Its orientation was corrected during derivation; early component-only status was later overtaken by other complete architectures. |
| Affine-chain parameterization | Adjacent affine frames reduced coefficient work. Complete synthesis, delayed frame phases, compact dynamic quotient creation, and streamed phase contractions passed. Width fell Q985 to Q955 at unchanged 951,995 T. | One-qubit lifetime coloring was not integrated; a dynamic remainder-step alternative increased controlled arithmetic. |
| Binary swap index | Six qubits replaced a bank of adjacent-swap flags; forced-swap and full native checks passed. | The approximate chart supports at most one swap; additional pivots contribute verifier errors. |
| Filled completion frames | Filling zero suffixes with the last active frame eliminated a boundary reciprocal. Full Q947 decoder passed with eight value/one phase errors. | A Q942 scan-reuse launch failed to import because its shell omitted the environment. Its empty JSON is not verification evidence. |
| Canonical and seedless affine synthesis | Removed a degree counter and endpoint field; completion ghosts and inverse phases were checked. An independently verified seedless implementation crossed milestone 11. | Equal active heads or accidental identity factors remain unsupported in the canonical chart. These approximations were charged rather than assumed absent. |
| Amplitude-polynomial preparation | Convert to `C=A/B' mod B`, absorbing fixed field factors into classical masks. The point scan drops from **80 to 11 T**, totaling **45,045 T** over 4,095 points. | Conversion and reversal are substantial and must remain in the total budget; a scan-only count is not an improvement. |
| Deferred converter measurements and inverse-mask kernels | Holding X-basis outcomes across the scan avoided retaining quantum intermediates. A degree-35 converter pair saved 12 Q and 3,541 T. Tested inverse-mask, 31-T trace, and 30-T joint/compact trace kernels passed phase/scratch checks. | Some compact trace variants remained standalone components; passing them does not mean every variant entered the final circuit. |
| Raw continued-fraction architecture | Removed scalar normalization and node swapping; degree-one convolution shares coefficient products. Complete small and full native oracles passed; first complete score above 13. | Rare extension poles persisted. Algebra-only supported-chart checks were not treated as all-input exactness. |
| Seedless raw endpoint and low-temporary synthesis | Reused source slots and outer polynomial relations; full circuits and fixed failure regressions passed, with preserved checkpoints around **13.4193** and **13.4855**. | Lower width sometimes increased gates; sample counts varied, so common failure cases were compared under two fixed measurement seeds. |
| Three-field terminal jet | `(p,q,d)` replaces four fields; tail `N=[d/q,p], B=[q]` makes the synthesized denominator monic. Removed both normalization passes: **557,645 -> 547,147 T** at Q896. | The physical zero-q chart required separate reversible handling. The independent run had five value/two phase errors; later destructive version superseded it. |
| Over-aggressive completion removal | Corrected implementation checked the preceding polynomial's leading ratio and passed padding/regression checks. | Initial counter-free chart mistook a genuine monomial gcd for a fake marker: 14 value/five phase errors and **765.5B** cost. The failing prototype was not selected. |
| Destructive quotient recovery | After recovering `head`, X-measure the redundant remainder `num=head*den`, reuse its 24 wires, and repay phase bilinearly. Quotient **664 -> 348 T**, private scratch **20 -> 16 Q**; complete circuit saved **12,572 T and four Q**. | An explicit 42-T correction for two absent proper-outer coefficients was necessary. The final tested circuit includes it. |

Final static DQI breakdown, which sums to the complete **534,575 T**:

| Stage | Toffolis |
| --- | ---: |
| Forward raw extraction | 114,240 |
| Forward rational synthesis | 34,750 |
| Forward amplitude conversion | 80,412 |
| Phase scan | 45,045 |
| Reverse amplitude conversion | 82,662 |
| Reverse rational synthesis | 45,376 |
| Reverse raw extraction | 132,090 |

Primary evidence: [final build](../initial_experiments/6/agent14/results/raw_cf_destructive_n70_m4095_build.json), [default run](../initial_experiments/6/agent14/results/raw_cf_destructive_n70_verify_default.json), [independent run](../initial_experiments/6/agent14/results/raw_cf_destructive_n70_verify_independent.json), and [destructive quotient source](../initial_experiments/6/agent14/destructive_quotient.py). The eight saved common failure cases matched the preceding implementations under two measurement seeds.

## Agent14: inverse-square-root experiments

| Experiment | Outcome |
| --- | --- |
| Adjacent Booth-product pairing | Three shared GF(2) products replace four elementary products while preserving the integer carries and rounding rules. Full verified Q89/T1727 circuit cost **175,049.68**. |
| Pairing at the rounded carry boundary | Further reduction to Q89/T1709, full cost **173,447.64**, unchanged approximation and 34 rotations. |
| Reconstruct curvature from live coefficients | An exact decoder using C0 low 21 bits, eight address bits, and slope high four bits shares 52 quadratic products plus an exception correction. Curvature can be released during linear multiplication: **Q85/T1763**, full cost **170,242.25**. This was agent14's final inverse. |
| Initial curvature replay schedule | Standalone decoder passed, but complete preflight found dirty scratch because live delayed C0 bits were borrowed as zeros. Corrected live-input decoder preserved the original timing. |
| Narrower coefficients | Curvature-7 fit passed the numerical domain with 324 leaves, but best tested complete schedule Q88/T1863 cost about **185,533**. Narrower slope/constant families failed numerical accuracy. |
| Exact coefficient-bank compression | The 47 coefficient bits, constant, and nine address bits had full binary linear rank 57; simple affine compression saved nothing. |
| Coefficient-bit embedding | Two curvature bits could be embedded with 324 leaves, enlarging lookup too much. Tested three-/four-bit constraints failed accuracy. |
| Rounded-square/fused correction | Replaced an unsuccessful 531-AND truth-table square with native **77/83-T** square and **93/99-T** full-coordinate components. Exhaustive coordinate, phase, and scratch checks passed. Numerical fused rounding with an extra fractional bit passed, but the complete prototype Q90/T2703 lost. Extra-bit-free local refits left infeasible leaves. |
| Cubic Horner fit | A 167-leaf fit passed numerical checks, but complete Q92/T1982 implementation lost. Narrower tested cubic families failed tolerance. |
| Radix-eight multiplication | Initial circuit failed: magnitude-three products require integer carries; direct XOR combination was invalid. No replacement selected. |

All three accepted Booth/curvature versions have full **200,540,928-shot** evidence. Agent20's EXP225a costs less than all of them, so the integration substitutes that whole verified component. No untested attempt was made to transplant agent14's arithmetic into EXP225a. Evidence: [final agent14 inverse run](../initial_experiments/6/agent14/results/pair_curve_live_q85_t1763_verify_default.json).

## Agent14: chemistry experiments

Agent14 retained the round-five R1/B719/C350 eight-bit structured physical factorization: Q454, T31,608, effective normalization **48.56733747**, cost **696,942,846.91**. Its compiler gain over Q478 belonged to round five, not new agent14 chemistry work.

| Experiment | Outcome |
| --- | --- |
| Direct B671 integer descent | Exact-word error **23.943463 -> 21.127105 mHa**; still above 20. |
| Full-angle Hessian/Newton proposals | 1,652 accepted row moves reduced the endpoint to **21.105792 mHa**; normalization about 56.31, also worse. |
| Lattice/LLL proposals | Reached **21.105777 mHa**; no feasible replacement. |
| Global QUBO, atom exchanges, shared-tree neighbors | Global proposals, 15 rounds of eight-atom exchanges, and 88 native-word tree neighbors did not improve the final endpoint. |
| B671 unrestricted coefficient-rank capacity | Full numerical PSD fit at the final fixed words remained **20.034993 mHa**. Increasing rank alone did not rescue those words; this is numerical evidence for that dictionary, not a bound on every B671 dictionary. |
| Reduce B719 coefficient rank | C340 passed tensor accuracy at about **19.910161 mHa**. Forty L1 Jacobi gauge sweeps only reduced normalization to **56.1249**, losing the small preparation saving. |
| Robust continuous and exact-forward STE fits | B511 continuous/refit **48.949/132.582 mHa**; straight-through continuation refit **166.686 mHa**. B671 STE and projected STE retained the original best but ended around **116.387/102.230 mHa**. |
| Quantization-aware ADMM and new geometries | R2/B383/C200 reached **19.379 mHa continuous** but about **156 mHa** after tested quantization continuation. R4/B255/C90 remained about **69 mHa** after continuation. B671 ADMM did not improve its integer initializer. |
| Disjoint dictionary splits | Six unrestricted block-covariance screens failed at **281–1,274 mHa**. |
| Overlapping dictionary splits | Twelve fixed dictionaries with shared atoms failed unrestricted span screens at **189–868 mHa**. |
| Exact factor-column partitions | Replicating U719 and splitting columns preserved modeled T while increasing Q; no gain. |
| One small dictionary per factor | Spectral truncation B7/B15 left about **4.771/2.985 Ha** error; retaining B76 reproduced the float incumbent and checked the implementation. |

None supplied a better native chemistry component. Evidence includes [B671 capacity](../initial_experiments/6/agent14/results/rank_capacity_671.json), [C340 gauge](../initial_experiments/6/agent14/results/rank_gauge_C340.json), [exchange endpoint](../initial_experiments/6/agent14/results/exchange_B671_final.json), and [tree endpoint](../initial_experiments/6/agent14/results/tree_B671_final.json).

## Agent20: successful experiments and adopted components

Agent20 began from round five's integrated **9.8110596433**, not agent14's final package. Its [final report](../initial_experiments/6/agent20/REPORT.md) supersedes old current-status text in [RESULTS.md](../initial_experiments/6/agent20/RESULTS.md) and [CONTINUATION.md](../initial_experiments/6/agent20/CONTINUATION.md).

| Experiment | Verified result and disposition |
| --- | --- |
| Local Booth/Karatsuba row rewrite | Exact arithmetic rewrite, unchanged approximation. Q90/T1693 cost **173,956.46**; Q89/T1709 cost **173,447.67**, each full default/independent PASS. Useful early gain, superseded by EXP225a. |
| Retained scalar inverses | First port saved gates but produced **50,501 phase errors/102,400**: source/reciprocal swaps invalidated retained norm intermediates. Correct source-preserving version passed small and both full runs at Q1004/T1,162,746, default cost **1.39613e12**. First milestone 10 followed. |
| 28-AND cube | Exhaustive arbitrary-target, phase, scratch checks passed independently, as in agent14's separate work. Generic C-fraction extraction improved, but this did not solve its zero charts. |
| Repaired Wronskian/Sylvester H endpoint | In-place bilinear division, padding-aware first-gap repair, and binary gap tag reduced point scan **80 -> 46 T**. Full Q1004/T1,123,942 decoder cost **1.3557830991e12**, default zero value/one phase error; independent zero/zero. Combined score **10.1289652741**. |
| Public Ahmed Berat Özer merged J-fraction decoder | Shared one reciprocal across two fraction levels; Q912 decoder passed two full runs, cost **1.0712306623e12**. With then-selected Phuc EXP222a inverse and local chemistry: **11.7473341475**. Saved original source/license/notice under [public/ahmed_v3](../initial_experiments/6/agent20/public/ahmed_v3/PORTFOLIO_NOTICE.md). |
| Public Phuc Thai EXP225a inverse and EXP119a decoder | Saved release bytes and provenance; full default/independent checks passed. EXP225a Q82/T1377/44 rotations, cost **138,095.0701**. EXP119a Q898, cost **1.0411801746e12**, zero observed errors. With unpaired chemistry: **11.933925**. Only the inverse is selected in this integration. |
| Paired chemistry ANGLE lookup | Adapted Phuc's fouraddr2 bilinear-combine idea to the existing compressed bank: `(a0+d2)(a1+d1)=a0*d1+a1*d2+d1*d2+a0*a1`. Two terms share a Toffoli; known products are compiled into the base table. Physical words/NPZ unchanged. **878 T saved**, 31,608 -> **30,730**, Q454 unchanged. Both full chemistry runs and all-address arbitrary-target/phase tests passed. Best packaged agent20 score **12.0465153528**. |
| Recovered gap flags | Derives flags from fraction heads and reverse polynomial numerators instead of storing a persistent flag bank. Q886/T1,100,012, full default/independent zero value/phase errors, costs **1.0233491626e12 / 1.0233492707e12**. Verified but not packaged there; loses to agent14's DQI. |

Recovered flags are explicitly approximate: the local zero/nonzero-head rule does not support all higher-order or consecutive ambiguous zero charts. Guard HMR clears exceptional residual bits. Forced unsupported charts fail; sampled zero-error runs do not establish an all-input inverse. The initial adapter's missing `km` binding was fixed before the successful frozen build.

Primary integration evidence: [EXP225a default](../initial_experiments/6/agent20/results/public_phuc_invsqrt_default_verify.json), [EXP225a independent](../initial_experiments/6/agent20/results/public_phuc_invsqrt_independent_verify.json), [paired chemistry default](../initial_experiments/6/agent20/results/paired_chemistry_default_verify.json), [paired chemistry independent](../initial_experiments/6/agent20/results/paired_chemistry_independent_verify.json), and [paired compiler](../initial_experiments/6/agent20/paired_chemistry_angles.py). Public sources were already saved and verified during agent20's experiments; this integration did not fetch or substitute a newer release.

## Agent20: DQI directions that did not become better selected components

| Direction | Result and limitation |
| --- | --- |
| Retention limits and fixed schedules | Limit 1004 pushed width to Q1016 and lost. Reoptimized one-/two-/three-slot schedule saved only a projected 340 T; no new full build. |
| Serial C fractions | Cube specialization reduced extraction **503,655 -> 466,583 T per direction**, complete static projection **1,450,238 -> 1,376,094 T**. Zero-head repair remained missing. |
| Compact/block/general Jacobi fractions | Two-pair sharing and larger affine blocks passed native components. Generic block-2 Q915/T1,104,604; block-3 Q987/T1,037,982; block-4 Q1035/T1,027,914. Natural-law zero-chart failures erased projected gains. |
| Leading-zero/variable-shift repair | Quantum-selected x²/x³ series inverse passed exact reciprocal/phase/involution checks. Repair reduced a tested model's wrong vectors, but missing later pivots remained. |
| Shared ordinary/gap synthesis | Small full oracle passed 4,096 sampled shots. Full profile Q942/T1,311,634 plus remaining exceptional-path penalty was uncompetitive. |
| Even/odd series reciprocal | N70 primitive **30,158 -> 27,762 T**; only 616 T saved in selected schedule, so no new full candidate. |
| Polynomial modular inverse by fraction reversal | Monic prototype Q868/T114,138; bilinear trimmed version Q892/T73,542; variable-prefix Q933/T78,882 forward. Native checks passed through degree 35, but interior gaps and complete routing/scan costs remained. |
| Projective modular inverse and in-place field reciprocal | Degree-35 projective Q939/T71,782 forward, native arbitrary-flag paired checks passed. In-place field reciprocal 66 T, controlled 91 T, with norm transport after source swap. Components, not a better full decoder. |
| Prefix counters, high-copy, compact head frames, lookahead, mixed quotients | Reusable components and small full model-phase oracles passed. Routed mixed projection Q943/T967,536 lost after a **658/100,000 wrong-vector** census; small oracle truth errors on unsupported charts were retained. |
| Triple-gap head/frame/block | Canonical head Q145/T2,498 per direction and integrated small blocks passed. N70 extraction alone **425,579 T per direction**, before synthesis/scan, was too expensive. |
| Generalized degree-1/2 quotient reconstruction | Native adjacent-flag, skip-priority, arbitrary paired-input tests passed. N35 v1 Q989/T52,454 per direction; improved v3 roughly Q954/T48,978. Reconstruction-only counts do not establish complete score gains. |
| Reversed-quotient cofactor/lookahead prefix | Components passed. N70 b6 block Q1000/T102,376 lost against the older generic block near 72k T. The report's `serial_matches: 2` compares different head encodings and does not independently diagnose tail failure. |
| Copied inverse/core profiles | Correct nested reciprocal-copy profiles lost. Earlier outer-only profiles did not exercise the copied inverse and cannot support a claimed gain. |
| Alternative phase endpoints and polynomial constructions | H, continuant, quotient, and polynomial components provided useful algebra/native evidence; no better full extraction/conversion/scan composition emerged beyond the recorded H/public/recovered candidates. |
| General multiplier/Forney decomposition searches | 59,395,257 tensor transformations in one bounded run retained rank 42; related quadratic searches found no cheaper selected joint decomposition. No mathematical optimality claim follows. |

Detailed component/version evidence is indexed in [head/quotient research manifest](../initial_experiments/6/agent20/results/head_quotient_research_manifest.json) and [generalized-quotient manifest](../initial_experiments/6/agent20/results/generalized_quotient_research_manifest.json). These are historical snapshots. Older recorder scripts assume older best hashes and should not be rerun unchanged against today's selection.

## Agent20: chemistry and inverse directions that lost

| Experiment family | Outcome |
| --- | --- |
| B719 lower-rank/gauge screens | C340 tensor accuracy was possible, but normalization about **57.24** (PCA gauge about 103) erased resource savings. |
| Quantization-in-loop, projected-gap and robustness controls | Initial unstable and smaller-step Adam/normalization-constrained runs did not improve changed integer words. Temporary fixed-word error improvement to 19.9375 mHa did not survive useful word changes. |
| B703 pruning/alternating repair | Integer descent and five-cycle alternation progressed through **20.414969/20.144103 mHa**. Later accuracy recovery still increased normalization/cost; no selected improvement. |
| Optimal-column partitions | Tested R2/B351 and R2/B479 dictionaries refit around **465/237 mHa**. |
| Fresh signed R4/B240/C88 fit | Continuous **63.5026**, rounded **173.9854**, exact-word PSD-refit **132.3009 mHa**, normalization about 78.70. Invalid. |
| Sparse positive-pair and density/Schur bootstraps | Sparse fixed-support floor **247.30 mHa**. Density target truncation 7.452 mHa looked promising, but decoded dictionaries refit to **591.083 mHa**, normalization about 140.35. |
| Conditioning-aware growth B479 -> B511 | Best **22.1378 mHa**, normalization 96.449; lowest normalization 71.453 occurred at **26.7321 mHa**. Invalid. |
| PSD generator factor via spectral/polar ADMM | Constructed 352 PSD generators at **7.36114 mHa** and intrinsic normalization about 47.91043. A valid matrix factor did not supply a compact integer angle dictionary. |
| Positive R4/B240/C88 dictionary fitting | Continuous improved through **38.083 -> 35.950 mHa**; integer refinement reached **56.729335 mHa**. Its saved union-span floor **34.541746 mHa** excludes weight-only rescue of those words. |
| Fixed-factor, weighted, global-atom, and longer-inner ADMM | Did not improve the saved integer initializer; terminal negative controls. |
| Positive R11/B64/C32 | **174.30 mHa continuous**, about **267.22 mHa exact-word NNLS**; fixed-dictionary span floor **158.634857 mHa**. |
| Signed continuation from positive R4/B240 | 40,000 steps reached **25.458428 mHa continuous / 136.921887 mHa exact-word refit**. Invalid; no selection change. |
| Rounded-square inverse | Exhaustive numerical model passed; complete Q122/T2063/35-rotation circuit passed 84,992 development shots with restored inputs/scratch but cost more. No full acceptance run warranted. |
| Fused-square, affine/static-top, split-bank and direct signed-row inverse variants | Numerical or native components worked, but full lookup, arithmetic, exponent routing, rotation synthesis, and lifetime costs lost. EXP225a remained the selected winner. |

The [signed-continuation result](../initial_experiments/6/agent20/results/signed_tensor_R4B240_from_positive.json) is a failed late experiment, not an improvement awaiting adoption. Necessary trace/rank constraints, fixed-dictionary spans, and finite fitting failures do not exclude every dictionary or factorization of a shape.

## Corrections, milestones, and remaining leads

The retained full-word field inverse costs **66 T**, Q45 including 12 source bits; early-release full-word version costs **102 T**, Q33 including source. **34 T is a specialized trace-inverse kernel**, not a full inverse. [Agent20's emitted-count correction](../initial_experiments/6/agent20/results/inverse_primitive_count_correction.json) supersedes informal 34/66 full-inverse estimates. Agent14's 30-T inverse phase repayment is likewise a specialized phase kernel, not a standalone complete reciprocal.

The chemistry normalization lower bound **46.78617959** holds under the fixed native definitions and physical domain. Relative to 48.56733747, normalization-only headroom is small. It does not bound circuit resources, other representations, or the combined score. Both agents' negative chemistry results repeatedly show that lower continuous residual, fewer rows, or lower rank alone is insufficient: decoded accuracy, normalization, and the full resource product must improve together.

Agent14 saved first crossings **10.2683683, 11.0283013, 12.0781345, 13.2024700, 14.0643040**, with hashes checked in [milestone audit](../initial_experiments/6/agent14/results/all_milestones_final_audit.json). Agent20 saved thresholds 10, 11 and 12 at **10.0304358, 11.7473341 and 12.0465154**. Their original archives and intermediate checkpoints remain intact. Integration establishes a score above 15 in its own directory; it does not rewrite historical first crossings.

Agent20's proposed same-target quadratic factoring for EXP225a remains **unimplemented and untested**. Factoring groups such as `CCX(a,s,t); CCX(b,s,t)` by reversible source parity changes is an exact-algebra lead, not a demonstrated additional reduction. Further decoder combinations likewise require a complete implementation, phase/scratch checks, and native acceptance before selection. The current delivered result includes only the already verified improvements.

This report consolidates the documented experiment families in both final reports and chronological logs. The linked source reports retain per-version measurements, failed logs, source hashes, and detailed reproduction commands. The [research evidence inventory](../initial_experiments/6/integration/results/research_inventory.json) lists the agents' saved JSON result artifacts and source-document hashes, including results not selected for packaging. The corresponding documentation copy is [optimization_iteration_6.md](optimization_iteration_6.md).
