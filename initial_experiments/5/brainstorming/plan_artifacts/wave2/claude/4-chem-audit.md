# Round 5 / chem-audit — bank the margin, then test chemistry's limiting assumptions

You own `initial_experiments/5/chem-audit`.

**Starting facts:**
- Incumbent: `initial_experiments/4/codex/candidates/norm_cycle2_signed`. Shape R1/B719/C350, 8-bit, J7, chi10.
- Resources: T=30,504 (inner prep 17,333, ANGLE 10,604, Givens 2,250, outer 317); Q = 141 persistent + 337 temporary.
- Physics: Λ=303.53, E_gap=3.915157 (floor 3.914898), λ_eff=48.594; error_2=19.9953 mHa against a 20 mHa cap.

Round-5 findings:
- `chem-qfit/discrete.py` takes the incumbent's words to 19.5164 mHa.
- Smaller passing fits cost more: B711 at 728.4M, B700/C360 at 825.3M.
- Fresh float basins at 13–15 mHa became 87–171 mHa after rounding and refit.
- `chem-blocks` multiblock fits never constrained Λ and came back with Λ = 985–1,838.
- Host-exact resource frontiers exist in `chem-blocks` and `chem-qfit`.

## Phase A — bank the margin (about a day)
Spend the freed ~0.48 mHa on λ_eff with constrained weight/offset refits (round-4 Codex method, AGY's fast offset solver as generalized in `chem-qfit/physical.py`). Native-verify. Replace the incumbent only if the cost improves.

## Phase B — assumption audit
1. Locate Codex's 46.78618 λ_eff bound in `initial_experiments/4/codex` and state exactly what it assumes. Does it bound every factorization under the fixed shift, or only the selected model?
2. Map the (error_2, Λ) frontier with fits from scratch, using Low et al.'s objective (arXiv:2502.15882 Eq. 41: residual + Λ penalty + gap term) at the verifier's fixed shift and exact metric. Cover the incumbent shape plus R1/B639, R2/B447–479, and one aggressive frontier shape.
3. Rounding robustness: train against random ±½-LSB angle perturbations (an expected-error objective, not straight-through), keeping weight norms bounded and Λ ≲ 320. Then run integer refinement with exact refits. Also try integer-native growth of a small shape.
4. At every checkpoint, track float error, rounded error, refit error and Λ.

## Success and scope
- Any native-passing candidate below 708.5M counts. Target ≤567M (1.25×), stretch ≤420M.
- Stop Phase B if the bound proves universal and no Λ-regularized float fit beats the incumbent's (20 mHa, Λ 303.5) point on the tested shapes. Document both findings.

## Shared constraints
- Write only inside your assigned directory. Do not read or search the other concurrent directories (`inv-width`, `dqi-dispatch`, `dqi-core`, `chem-audit`). Everything else may be read; `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Start from the round-4 Codex-selected artifacts. Local reference **9.3677221531**: inverse 182,956.5272, DQI 1,598,204,229,385, chemistry 708,547,388.9553. Preserve them until a replacement passes.
- The public leaderboard top is **18.9**, an 8.2× larger product of task ratios than ours, so better solutions exist. A failed approach is evidence about that design only; report it that way.
- Score with the actual objective: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Develop on sampled inputs and short runs, but accept only on full default native verification: DQI needs all 102,400 shots plus independent salts; chemistry needs the native ANGLE/oracle and physical checks. No hash or salt grinding.
- Four agents share 8 cores and one GPU. Respect your worker cap, and wrap full-default verifier runs in `flock /tmp/r6-fullverify.lock …` (the only permitted write outside your directory).
- Report a numeric first checkpoint within a few hours. Maintain `RESULTS.md`: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.