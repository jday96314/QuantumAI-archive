# Round 5 / chem-qfit — quantization-aware factorization compression

You own `initial_experiments/5/chem-qfit`. Focus on recovering **a smaller physically valid eight-bit chemistry factorization** by making integer-angle fitting materially better. Do not merely continue float-only L-BFGS or optimize normalization of the incumbent.

**Starting facts:** Round-4 Codex's passing chemistry submission is `R=1, B=719, C=350, brot=8, J=7, chi=10`, signed-angle lookup, `T=30,504`, `Q=478`, `lambda_eff=48.59417497035`, cost **708,547,388.9553**, one-body error **11.986315 mHa**, two-body **19.9953045 mHa**. The latter has almost no margin. The actual bottleneck is inner state preparation (17,333 Toffolis and 337 temporary qubits). Claude found promising float fits for some smaller two-block shapes, but the attempted local quantizer turned ~13 mHa into ~25 mHa; its control also failed to reproduce an existing good quantized basin.

## Work plan

1. **Calibrate before innovating.** Recover the exact accepted round-4 integer angle words, tree, weights, and offsets, and reproduce their physical metrics using the same reference tensors, float32 conventions, symmetry shifts, angle encoding and signed gauges as the native verifier. Demonstrate that decoding/re-encoding the incumbent leaves the solution valid. Reproduce the current integer-angle refit (round 3/4) as a control; characterize why Claude's quantizer lost accuracy on a known accepted basin.

2. Develop at least two **meaningfully different discrete optimizers**: e.g. Givens-tree error-feedback rounding with residual refits, progressive vector freezing with retraining, and small-block/coupled integer-angle moves with analytical weight/offset refits. Optimize the true two-body tensor error and one-body residual, not mean angle error or an unquantized proxy. Try efficient batched evaluation and use AGY's fast offset solver only after confirming its exact semantics.

3. Use the **actual host cost model**, including chosen QROAM exponent `k2`, inner temporary peak, full selected signed-angle oracle costs, and normalization, to choose a handful of **small-table shapes**. Include an 8-bit single-block `B~511–700` corridor and the previously investigated two-block `B=479/447` corridor with realistic `C`; test `R2/B351/C175` only if a continuous fit has sufficient margin. Distinguish `B per ensemble` from total `R*B`. Do not assume a shape scores well from row count alone.

4. Run matched comparisons of continuous fit, direct rounded fit, and **refitted exact integer-word fit** for each candidate. First establish whether a shape can reach <=20 mHa with enough slack; then build complete ANGLE circuits and verify the host result. The physical checks include 20 mHa on both one-/two-body residuals, the strict Lambda/gap floors, XOR-difference oracle interface, <=128 qubits per oracle, and correct reversal.

## Success criteria and scope limit

- Preferred complete win: chemistry cost <=**~567M** (>=1.25x gain), stretch <=**~506M** (>=1.40x), with native passing verification. Preserve any genuinely improved quantizer and numerical diagnostics even if no submission passes.
- A float-only pass or a passing incumbent with worse cost does not count. If smaller shapes repeatedly fail, compare on an identical float basin and classify failure (shape capacity vs numerical optimizer vs quantizer) instead of declaring a general floor.

## Shared constraints

- Write only inside your assigned experiment directory. **Do not read, search, or copy any other round-5 agent directory.** Prior rounds, shared baseline assets, and `from_host/` are read-only. Do not modify the verifier or official scoring rules. Discover actual repository paths rather than assuming links in summary reports are present verbatim.
- Start with the **round-4 Codex-selected**, not round-3 or another agent's inherited, artifacts. The local combined reference is **9.3677221531**: inverse cost **182,956.5272**, DQI **1,598,204,229,385**, chemistry **708,547,388.9553**. Preserve those candidates until a replacement passes.
- Use the **actual objective** rather than a gate-only proxy. For DQI: `R(epsilon_plus) * Q**2 * mean_Toffolis` using the host's Clopper–Pearson definition; for chemistry: `lambda_eff * T_total * Q_total`; for inverse: `Q * (mean_Toffolis + total_T/4)`.
- Maintain `RESULTS.md` with tested hypotheses, implementation changes, checks, exact/verifier-measured vs estimated numbers, failed candidates, and clear reproduction commands. Deliver any winning circuits/NPZ, builders, and native verifier JSONs. A faster component or a floating-point fit is **not** a verified task-cost win.
- For final DQI acceptance, use all **102,400** shots and independent validation salts; for chemistry, complete the native ANGLE/oracle and physical validation (not just the heuristic resource model). Do not exploit canonical-hash/sample variation, alter verifier confidence, or claim a win from shortened tests. Keep independent controls and attribute reused Apache-2.0 public code.