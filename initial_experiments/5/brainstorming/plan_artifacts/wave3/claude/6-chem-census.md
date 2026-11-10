# Round 5 / chem-census — chemistry routes beyond the incumbent basin

You own `initial_experiments/5/chem-census`. Goal **≥1.25×** (≤ ~566M).

**Facts:**
- λ_eff ≥ 46.786 is proven for every factorization at the fixed shift (`initial_experiments/5/chem-audit/RESULTS.md`).
- Banked incumbent: 708.156M; λ_eff ≈ 48.567; T 30,504 (inner prep 17,333 ≈ R·C·⌈(B+1)/16⌉ + small terms, ANGLE 10,604, Givens 2,250, outer 317); Q 478 = 141 + 337. The host's k2 rule pins Q near 460–480, so **T is the lever**.
- Evidence from `chem-qfit`, `chem-blocks` and `chem-audit`: the incumbent's integer-native basin sits at 19.5–20 mHa; Λ-regularized fits from scratch reach only 28–47 mHa before rounding; freshly quantized float basins lose 74–156 mHa; multiblock fits failed.

**Phase 1:** propose ≥5 distinct routes (not repeated shape sweeps), each with an objective estimate and its biggest uncertainty. Consider:
- continuation from the incumbent integer basin to smaller tables, crossing B+1 multiples of 16 (each step cuts ~560 T, so λ_eff may rise <1.8% per step);
- ANGLE oracle designs that exploit dictionary structure;
- reallocating the unused one-body error budget (12 of 20 mHa);
- the host model's discrete thresholds mapped onto continuation paths;
- integer-native growth of a smaller shape.

**Phase 2:** test the top two routes' biggest uncertainty with native physical evaluation. **Stopping** as in the other censuses. GPU plus 2 CPU workers; coordinate GPU use with chem-audit if it is still running.

## Shared constraints
- Write only inside your assigned directory. Reading other directories is allowed but read-only; treat other agents' in-progress files as unstable. `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Current local reference: combined 9.3677 (inverse 182,956.5272; DQI 1,598,204,229,385; chemistry 708,156,073.7876 banked in chem-audit). Preserve selections until a replacement passes.
- The public leaderboard top is 18.9 (≈8.2× our product of ratios). A failed approach is evidence about that design only.
- Objectives: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Accept only on full default native verification (DQI: all 102,400 shots plus an independent salt; chemistry: native ANGLE/oracle and physical checks; inverse: all 200,540,928 shots). No hash or salt grinding.
- The machine is shared (8 cores, 1 GPU): respect your CPU cap and wrap full verifier runs in `flock /tmp/r6-fullverify.lock`.
- Report a numeric first checkpoint within a few hours. Keep RESULTS.md current: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.