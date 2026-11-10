# Round 5 — promising directions outside the top five

**Decision criterion.** These ideas have plausible technical merit, but I would *not* allocate one of the initial three round-4-sized Codex efforts to them. In most cases, either (a) their attainable standalone task gain appears modest, (b) they depend on progress from a selected top-five experiment, or (c) the existing evidence is too preliminary to justify a full agent run. None of the exclusions below is a proof of infeasibility.

**Baseline:** selected round-4 local combined score **9.3677221531**; `invsqrt` 182,956.5272; DQI 1.598204229385e12; chemistry 708,547,388.9553. DQI minimizes `R(eps_plus) * Q**2 * T`, chemistry minimizes `lambda_eff * T * Q`, and inverse square root minimizes `Q * (T + synthesized_T/4)`. These are *scored* resources, not simulator runtime.

## Secondary ideas worth retaining

| Direction | Key experiment / falsifiable checkpoint | Why below the top five |
|---|---|---|
| **1. Content-aware ANGLE lookup synthesis** (`chem-qrom`) | On the **exact passing round-4 angle words**, compile a content-pruned trie / shared multi-output Boolean circuit. Exploit repeated subtrees, identical bitplanes, cheap Clifford frames, row reordering *with corresponding weight permutation*, and signed-angle gauges. Verify forward/reverse XOR-difference semantics and the <=128-qubit cap. Compare *measured* ANGLE Toffolis to 10,604 for two passes. | ANGLE is ~35% of the current chemistry Toffolis and not the scored-width bottleneck. A 50% ANGLE-only saving gives only ~1.21x chemistry. Important as a follow-on after shape compression, where lookup's fraction could rise. |
| **2. Dedicated Givens-tree and orbital-order optimization** (`chem-tree`) | On a fixed accepted physical dictionary, test tree changes and common orbital permutations. Optimize the integer-angle **tensor residual**, angle-word correlations, and compiled lookup complexity jointly. Require a valid recovered NPZ and tree; compare against the accepted tree with the same precision. | Closely coupled to `chem-qfit` (rank 2). Standalone numerical/encoding improvements may not move the inner-QROAM width. Avoid two agents independently rebuilding the same quantizer. |
| **3. Finish the Wronskian / half-GCD rational-reconstruction decoder** (`dqi-wronsk`) | Build a complete reversible pipeline around Codex's component that reduced one padding adapter from 28,069 to 7,140 Toffolis per direction. Include the evaluator, root scan, inverse, liveness and exceptional-case behavior; then calculate full `R*Q²*T`. | The saving is **component-level**, not yet end-to-end. Rational reconstruction can demand large intermediate buffers at degree 70. Escalate if `dqi-phase` identifies a compatible low-width dataflow. |
| **4. Bounded coefficient windows and reversible history reconstruction** (`dqi-windows`) | Profile coefficient support under >=1M law-sampled syndromes. For each step, compute the probability of exceeding a truncated window, and verify a shared-path update where displaced discrepancies/history are reconstructed from an invariant rather than saved. Measure gates and qubit peaks. | Published precedent is strong, but savings within a single existing compact algorithm may be incremental. Components belong naturally in `dqi-hybrid` after it establishes a viable storage schedule. |
| **5. Linear phase extraction from classical DQI masks** (`dqi-mask`) | Adapt Gu's syndrome-linear parity correction to **per-shot classical** masks. Derive a syndrome-controlled Clifford phase that removes the need to scan certain positions without assuming the mask is a build-time constant. Price all classical control/Clifford overhead and verify arbitrary masks. | Primarily a subproblem of `dqi-phase`. The published 71-position skip is not directly reusable under the Kaggle mask interface. Maximum standalone Toffoli saving is limited unless it enables a larger scan redesign. |
| **6. Inverse-square-root width-first scheduling** (`inv-width`) | Profile the selected 90-qubit circuit's exact live-register peak. Test in-place normalization, borrowing coefficient registers as accumulators, measurement-erased carries, and low-ancilla adders under an exact integer model. Require all 3,133,452 inputs and 64 repetitions in final verification. | Plausibly ~1.1–1.3x inverse-task gain, but the prior four rounds already optimized workspace aggressively. A useful fallback if DQI or chemistry opportunities fail, not a first-choice 17+ driver. |
| **7. Adaptive-prefix piecewise-linear inverse-square-root approximation** (`inv-linear`) | Fit exact integer acceptance windows for nonuniform mantissa intervals and compare one short linear product plus extra QROM against the current two-product quadratic. Evaluate QROM leaf cost and erasure, worst-case synthesis margin, and width. | The public Gu Coulomb example has **much looser** relative accuracy than this task. At 2^-24 relative error, lines may require too many intervals; direct quadratic fitting has already been particularly successful. |
| **8. Chemistry fast atom generation / offset refits** (`chem-seed`) | Reproduce AGY's reported 4.4-second analytical offset fit and tensor-power seeding; use them to create many *distinct* starting basins for `chem-qfit` or `chem-blocks`. Track quantized reconstruction **and** lambda_eff, not only a float residual. | Enabling infrastructure, not a standalone quantum-resource breakthrough. AGY's reported residual-floor improvement (19.4272 -> 19.1568 mHa) did not yield a better packaged submission. |
| **9. Chemistry normalization beyond the incumbent parametrization** (`chem-norm`) | Recheck the scope of the reported `lambda_eff >= 46.78618` bound. Explore transformations not covered by that *fixed-shift* bound only if they respect the verifier's ignored/immutable symmetry-shift fields and strict `Lambda`, `E_gap` constraints. | If the reported bound applies, maximum normalization-only gain is ~3.9%. The verifier fixes symmetry shifts; editing submitted `H1_I` / `Hsym_lower` does not change those physically. |
| **10. DQI alternative error-budget allocations** (`dqi-budget`) | Vary number of decoding iterations, coefficient cutoffs, and exceptional-path capacity. Estimate expected error rate on >=1M law-sampled draws and price its **one-sided 2^-128 Clopper–Pearson** repetition multiplier. | Could give a few percent or enable a compressed decoder, but rare failures are already expensive under `R(eps_plus)`. A 380-error architecture is a poor base for deliberately adding more approximation. |

## If one of the initial agents finishes early

**Fastest self-contained backup** beyond the selected top five: run `chem-qrom` directly on the *unchanged* accepted ANGLE words. It does not depend on discovering a new valid physical fit, so it can produce a decisive measured cost comparison quickly. Give it a strict go/no-go checkpoint based on measured **full-oracle** Toffoli savings, including reversible cleanup and the 128-qubit cap.

**Next-lowest integration risk:** `inv-width` on the selected round-4 90-qubit circuit. It has a smaller expected score gain, but every successful candidate can be tested independently of the chemistry/DQI branches.

**Most speculative high-upside backup:** `dqi-wronsk` if the division-free scan work uncovers an economical way to share/recycle polynomial storage.

## Cross-cutting instructions for any later experiment

- Do not change or bypass `from_host/` verification, reduce final shot counts, alter confidence, or select candidates by lucky salts/circuit hashes.
- Label separately: exact mathematical result, float fit, *quantized* physical fit, circuit component count, complete native pass, and actual combined-score improvement.
- For chemistry, never claim a resource-only shape is feasible before checking the recovered quantized angle words, 20 mHa coefficient residuals, fixed shifts, `Lambda/E_gap` constraints, XOR-difference oracles, and host-selected QROAM exponent.
- For DQI, never assume a quantum-controlled exceptional path is gate-free on ordinary inputs, or that a published *fixed* mask is equivalent to Kaggle's per-shot classical masks.
- Preserve generated candidates, full verification JSON, scripts, liveness/width reports, command lines, and a concise result ledger.

## Sources and context

Local project summaries supplied for this planning: `optimization_iteration_3.md`, `optimization_iteration_4.md`, `REPORT(2).md`, `claude_thoughts.md`, `verifiers.md`, and the selected `norm_cycle2_signed_verify.json`.

Public technical reference implementations:
- Gu et al., DQI Reed–Solomon:
  https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md
- Gu et al., Fe₂S₂ content-pruned lookup:
  https://andigu.github.io/projects/uncompute/chemistry-fe2s2.html
- Gu et al., Coulomb approximation:
  https://andigu.github.io/projects/uncompute/coulomb-potential.html