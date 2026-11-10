# Round 5 / chem-blocks — resource-aware multiblock factorization

You own `5/chem-blocks`. This is a **backup** exploratory architecture search, separate from the quantizer-improvement agent. Attack the chemistry **inner preparation table and qubit-width bottleneck** through structurally different factorizations, not merely deleting atoms from R1/B719/C350.

**Starting facts:** the accepted round-4 chemistry shape is R1/B719/C350, 8-bit, J7, chi10, `T=30,504`, `Q=478`, cost 708.55M. The host sets the inner QROAM exponent automatically using the angle-storage threshold. For example, `R4/B180/C88` has an attractive *conditional* 7-bit resource estimate because the lower precision can force `k2=4` rather than the 8-bit `k2=5`, reducing the inner difference-bank workspace by hundreds of qubits. No accurate 7-bit fit of that shape has been demonstrated, and earlier grouped fits failed badly after quantization. It is a candidate, not a prediction.

## Work plan

1. Build an **exact host-resource frontier** before expensive fitting: scan `(R,B,C,brot,J)` with actual `from_host` QROAM choice (`optimize_qubits=True`, not a custom selected `k2`), persistent+maximum temporary widths, analytically priced state preparation and Givens gates, and measured or explicitly bounded costs of the round-4 signed-gauge ANGLE oracles. Include 7/8/9-bit angles and discrete width cliffs near B=256/257. Identify the most promising valid-cost *shapes*, but label physical feasibility untested.

2. Derive actual symmetric generator matrices for the incumbent's two-body factorization and verify reconstruction with the **symmetric feature Gram metric** (do not substitute naive Euclidean PSD projection). Explore orthogonal mixing of generator columns, clustering them by approximate simultaneous diagonalizability, then **joint approximate diagonalization** for per-ensemble bases. Rotate associated offsets consistently. This is distinct from AGY's orbital-disjoint block split, which discarded inter-block Coulomb correlations.

3. Fit/relax one or two attractive shapes, considering `R2/B351/C175` at 8 bits and `R4/B180/C88` at 7 bits **only if their resource and continuous-fit screens justify them**. Permit overcomplete, nonorthogonal vectors when joint diagonalization is too restrictive. Use warm starts, constrained weight/offset refits, correct physical tensor metrics and a factorization-aware angle quantizer; never conclude 7-bit impossibility from one rounded fit.

4. Check the **actual quantized** one-/two-body thresholds (both <=20 mHa), strict Lambda and E_gap constraints, angle XOR-difference interface and 128-qubit per-oracle cap. Integrate and native-verify only promising candidates; report λ_eff, T, Q, total cost and all accuracy margins.

## Success criteria and scope limit

- Aim for **>=1.5x chemistry task improvement**, stretch near 2x, from a complete valid submission. A resource-only 2x projection is not a success.
- If the proposed R4/B180 seven-bit family cannot approach the required residual under reasonable numerical controls, pivot once to the best nearby **resource-favorable 8-bit** family, document results and stop rather than running unbounded float optimizations.

## Shared constraints

- Write only inside your assigned experiment directory. **Do not read, search, or copy any other round-5 agent directory.** Prior rounds, shared baseline assets, and `from_host/` are read-only. Do not modify the verifier or official scoring rules. Discover actual repository paths rather than assuming links in summary reports are present verbatim.
- Start with the **round-4 Codex-selected**, not round-3 or another agent's inherited, artifacts. The local combined reference is **9.3677221531**: inverse cost **182,956.5272**, DQI **1,598,204,229,385**, chemistry **708,547,388.9553**. Preserve those candidates until a replacement passes.
- Use the **actual objective** rather than a gate-only proxy. For DQI: `R(epsilon_plus) * Q**2 * mean_Toffolis` using the host's Clopper–Pearson definition; for chemistry: `lambda_eff * T_total * Q_total`; for inverse: `Q * (mean_Toffolis + total_T/4)`.
- Maintain `RESULTS.md` with tested hypotheses, implementation changes, checks, exact/verifier-measured vs estimated numbers, failed candidates, and clear reproduction commands. Deliver any winning circuits/NPZ, builders, and native verifier JSONs. A faster component or a floating-point fit is **not** a verified task-cost win.
- For final DQI acceptance, use all **102,400** shots and independent validation salts; for chemistry, complete the native ANGLE/oracle and physical validation (not just the heuristic resource model). Do not exploit canonical-hash/sample variation, alter verifier confidence, or claim a win from shortened tests. Keep independent controls and attribute reused Apache-2.0 public code.