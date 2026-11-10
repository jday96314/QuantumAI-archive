# Round 5 / dqi-phase — low-width division-free DQI phase scan

You own `initial_experiments/5/dqi-phase`. Replace the expensive **phase extraction** of the round-4 compact DQI decoder while preserving its compact storage. **Do not work on exception repair or reimplement the locator recurrence**; another agent owns that.

**Starting facts:** round-4 DQI has `Q=961`, `T≈1,019,866`, `R≈1.69685`. Its position kernel costs ~88 Toffolis over 4,095 field positions (~360k gates). Gu et al. publish a 23-Toffoli point body (11-root test + 12 pairing), but their 4,024-position scan and build-time parity simplification assume a **fixed mask**. Kaggle supplies per-shot **classical masks**, and Gu's derivative-inverse preparation uses large cofactor storage. A straight transplant can lose under `Q**2`.

## Work plan

1. Reproduce the current phase kernel's functional behavior and gate/width attribution. Read `from_host/verifiers/dqi/` and Gu's disclosed DQI architecture. Establish field basis/position orientation, the exact phase identity and mask dependence with an **independent classical oracle**; test randomized squarefree locators of degrees 0–35, exceptional cases and arbitrary classical masks.

2. Measure a small **near-term control**: can the published 80-Toffoli division-based point logic or a basis transformation beat the 88-Toffoli kernel at unchanged width? Retain it only if the full scan/circuit improves.

3. Main experiment: compute a reciprocal derivative polynomial `Psi` satisfying `Psi * sigma_prime = 1 mod sigma`, then either (a) scan via Gu's direct bilinear trace kernel, or (b) replace evaluator `omega` in place by `F = omega * Psi mod sigma`, so `F(root)` is the error magnitude. The `F` approach might reduce the pointwise nonlinear work to the root test plus Clifford parity/phase operations, **but this has not been demonstrated**. Derive and verify its exact cost with per-shot classical masks, including conditional phase, any 71-point skip, and all translation/correction overhead.

4. **Most important:** design and implement a reversible workspace/liveness schedule for building, retaining and erasing `Psi` or `F`. Include forward and reverse polynomial transforms, EEA/half-degree Euclid or alternative, all scratch, and complete restoration of 70 syndromes. Do not count a mathematically invertible map as a free reversible implementation. Profile actual live-qubit peak and source of every extra word.

5. Before full integration, produce a measured component budget and an **end-to-end** cost estimate under `R*Q**2*T`. If it beats the incumbent by >=15–20% with realistic overhead, integrate and native-verify with the current decoder. Otherwise document the least expensive working kernel and a precise obstruction.

## Success criteria and scope limit

- Aim for **>=1.25x full DQI cost improvement** on the existing failure distribution; meaningful component gains may be banked even if below that threshold.
- Keep width close to 961 if possible, but use actual break-even arithmetic—not an arbitrary +100-qubit prohibition. Never claim 11 or 23 Toffolis/point without all required extra preparation, inverse work, mask handling, and root gating.

**External reference:** https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md

## Shared constraints

- Write only inside your assigned experiment directory. **Do not read, search, or copy any other round-5 agent directory.** Prior rounds, shared baseline assets, and `from_host/` are read-only. Do not modify the verifier or official scoring rules. Discover actual repository paths rather than assuming links in summary reports are present verbatim.
- Start with the **round-4 Codex-selected**, not round-3 or another agent's inherited, artifacts. The local combined reference is **9.3677221531**: inverse cost **182,956.5272**, DQI **1,598,204,229,385**, chemistry **708,547,388.9553**. Preserve those candidates until a replacement passes.
- Use the **actual objective** rather than a gate-only proxy. For DQI: `R(epsilon_plus) * Q**2 * mean_Toffolis` using the host's Clopper–Pearson definition; for chemistry: `lambda_eff * T_total * Q_total`; for inverse: `Q * (mean_Toffolis + total_T/4)`.
- Maintain `RESULTS.md` with tested hypotheses, implementation changes, checks, exact/verifier-measured vs estimated numbers, failed candidates, and clear reproduction commands. Deliver any winning circuits/NPZ, builders, and native verifier JSONs. A faster component or a floating-point fit is **not** a verified task-cost win.
- For final DQI acceptance, use all **102,400** shots and independent validation salts; for chemistry, complete the native ANGLE/oracle and physical validation (not just the heuristic resource model). Do not exploit canonical-hash/sample variation, alter verifier confidence, or claim a win from shortened tests. Keep independent controls and attribute reused Apache-2.0 public code.