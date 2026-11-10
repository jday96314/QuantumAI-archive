# Round 5 / dqi-hybrid — compact Berlekamp–Massey / Euclidean hybrid

You own `initial_experiments/5/dqi-hybrid`. This is a **backup** research effort: find a decoder with the inexpensive shared-path locator updates of round-3/Gu Berlekamp–Massey and the much smaller live storage demonstrated in round 4.

**Starting facts:** Round 3 achieved Q=1,429, ~645,755 Toffolis, 8 phase errors/102,400, R≈1.223, DQI cost 1.613T. Round 4 achieved Q=961 but T≈1.020M, 380 phase errors, R≈1.697, DQI cost 1.598T. Gu's published shared-path BM construction is ~442k gates but Q=1,837 for its entire pipeline. The cheap components and low width have **not** yet been achieved together.

## Work plan

1. Identify where the round-3 BM and round-4 compact Euclidean decoder keep locator, template, cofactor, syndrome, evaluator, reversibility/history records and phase-kernel scratch. Generate comparable **live-register timelines** and forward/reverse Toffoli breakdowns. Start from the exact code that generated the selected submissions, not summary estimates.

2. Construct a **concrete representation** that reuses 840 syndrome-input qubits for BM and evaluator/cofactor state as soon as a syndrome is mathematically consumable, or a similarly tight reversible streaming schedule. Use Gu's shared arithmetic path, reversible branch swap, invariant-based recovery of discrepancy and bounded coefficient windows where genuinely compatible with the host's dynamic classical masks.

3. Prove each proposed overwriting/reconstruction operation is reversible on all tested inputs, including exceptional locator degrees, and confirm all syndrome bits return unchanged and scratch clears. An abbreviated liveness spreadsheet is not an implementation: emit a working component with native gate counts.

4. Implement the most favorable hybrid, with the **existing** phase scan initially to isolate the decoder gain. Validate under 1M law-sampled classical cases where practical and native default DQI shots before evaluating possible integration of the cheaper phase kernel.

## Success criteria and scope limit

- Target a complete circuit around **Q<=1,050, T<=600k, R<=1.23**; this would yield roughly a 2x DQI gain. Lesser complete gains are valuable.
- If no credible schedule below ~1,200 live qubits emerges after a full register-liveness/design pass, stop before rewriting the entire decoder and document the specific information-theoretic/reversibility obstruction (if proven) versus implementation difficulty (if not).
- Do not equate Gu's lower gates at 1,837 qubits with a competition win; price every candidate under **Q squared**.

**External reference:** https://github.com/andigu/uncompute/tree/main/tasks/dqi-reed-solomon

## Shared constraints

- Write only inside your assigned experiment directory. **Do not read, search, or copy any other round-5 agent directory.** Prior rounds, shared baseline assets, and `from_host/` are read-only. Do not modify the verifier or official scoring rules. Discover actual repository paths rather than assuming links in summary reports are present verbatim.
- Start with the **round-4 Codex-selected**, not round-3 or another agent's inherited, artifacts. The local combined reference is **9.3677221531**: inverse cost **182,956.5272**, DQI **1,598,204,229,385**, chemistry **708,547,388.9553**. Preserve those candidates until a replacement passes.
- Use the **actual objective** rather than a gate-only proxy. For DQI: `R(epsilon_plus) * Q**2 * mean_Toffolis` using the host's Clopper–Pearson definition; for chemistry: `lambda_eff * T_total * Q_total`; for inverse: `Q * (mean_Toffolis + total_T/4)`.
- Maintain `RESULTS.md` with tested hypotheses, implementation changes, checks, exact/verifier-measured vs estimated numbers, failed candidates, and clear reproduction commands. Deliver any winning circuits/NPZ, builders, and native verifier JSONs. A faster component or a floating-point fit is **not** a verified task-cost win.
- For final DQI acceptance, use all **102,400** shots and independent validation salts; for chemistry, complete the native ANGLE/oracle and physical validation (not just the heuristic resource model). Do not exploit canonical-hash/sample variation, alter verifier confidence, or claim a win from shortened tests. Keep independent controls and attribute reused Apache-2.0 public code.