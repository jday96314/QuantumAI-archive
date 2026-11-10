You own `initial_experiments/5/inv-census`. Goal **≥1.5×**: cost Q×(T + T_count/4) ≤ ~122k, e.g. Q 70 × 1,740, Q 60 × 2,030, or Q 80 × 1,525.

**Incumbent:** Q 90, mean T 1,782, 34 rotations (251 Toffoli-eq), cost 182,957.
- Stage costs: normalization 210, coefficient load + erase 304, Horner 1,046, output shift 217, affine restore 16.
- The 90-qubit peak is the final Horner product's first radix-4 row: 22 input, 5 exponent flags, 23 constant remainder, 17 slope, 4 curvature, 10 output headroom, plus carries.

**Already tried (read `initial_experiments/5/inv-width/RESULTS.md`; don't repeat without a new twist):** 423 incumbent schedules (best 88 Q at 183.2k); segmented phasing (the carry spans 6–9 bits); dirty-input borrowing (−12 Q for +1–2k T); three first-draft families at 448k–2M.

**Phase 1 (no building):** propose ≥6 distinct architectures, each with a peak register table, a T estimate, rotation count and synthesis margin, the objective, and its biggest uncertainty. Include:
- removing the output shift and its 10 headroom qubits (apply the exponent through the coefficients, the rotations or the normalization instead);
- different approximation forms with fewer coefficient bits or narrower live registers;
- in-place exponent handling;
- Gu et al.'s numerical-task techniques adapted to 2⁻²⁴ relative accuracy: power-of-two-aligned intervals, a shared local coordinate, narrow direct accumulation, measurement-based table erase with phase repair.

**Phase 2:** prototype the biggest uncertainty of the top two with exhaustive integer models and sampled native runs. **Phase 3:** a build plan for anything projecting ≥1.3×.
**Stopping:** never stop because a first draft loses to a four-round incumbent; drop a candidate only when its key uncertainty is resolved against it. CPU cap 3. Full verification is 200.5M shots, so develop with exhaustive integer models.

## Shared constraints
- Write only inside your assigned directory. Reading other directories is allowed but read-only; treat other agents' in-progress files as unstable. `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Current local reference: combined 9.3677 (inverse 182,956.5272; DQI 1,598,204,229,385; chemistry 708,156,073.7876 banked in chem-audit). Preserve selections until a replacement passes.
- The public leaderboard top is 18.9 (≈8.2× our product of ratios). A failed approach is evidence about that design only.
- Objectives: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Accept only on full default native verification (DQI: all 102,400 shots plus an independent salt; chemistry: native ANGLE/oracle and physical checks; inverse: all 200,540,928 shots). No hash or salt grinding.
- The machine is shared (8 cores, 1 GPU): respect your CPU cap and wrap full verifier runs in `flock /tmp/r6-fullverify.lock`.
- Report a numeric first checkpoint within a few hours. Keep RESULTS.md current: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.