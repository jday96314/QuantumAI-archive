# Round 5 / dqi-core — attribute and shrink the compact Euclid core

You own `initial_experiments/5/dqi-core`. This is research: you are **read-only** on the selected decoder and work on copies in your directory. A parallel agent owns the full decoder rewrite.

**Starting facts:**
- The compact core costs 649,802 T: 324,901 each way. The scan adds 370,152.
- Round-3 BM costs 293,503 (159,664 forward + 133,839 reverse).
- Gu et al.'s BM reverse is 61% of its forward, because it measurement-erases forward temporaries and recomputes only what the reverse needs. Our reverse equals our forward exactly.
- Compact state lengths are (k−1, k, 70−2k, 1) words. Blocks of three use a 6-word head; updates go serial after k=22. Dense triangular (tail-series) transforms appear to dominate.
- Reuse the replay and liveness tooling in `initial_experiments/5/dqi-hybrid` and `initial_experiments/5/dqi-phase/baseline_audit.py`.

## Work plan
1. Per-operation Toffoli attribution for both directions: series division and middle products, cofactor updates, six-moment transfer computation, inverses and scaling, flag/completion logic, padding.
2. List the forward temporaries that could be measurement-erased and recomputed selectively in reverse, and estimate the asymmetric-reverse saving.
3. Algorithmic options:
   - Replace per-step tail-series division with a recurrence on running remainder products.
   - Use larger half-gcd-style blocks early, where width is lower (trading head width against block size).
   - Find a decoder whose end state is (σ, ω) or (σ, F) within the same 852 qubits, enabling the 80-T Forney body or a cheaper scan.
4. Prototype the best option as a native component with exact counts and reversibility tests. Estimate the end-to-end cost under R·Q²·T assuming R≈1.2.

## Success and scope
- Deliver a design memo plus prototypes showing a ≥20% core reduction at unchanged peak width (≈1.15× DQI). Stretch: ≥30%, or a (σ, ω) end state.
- If the attribution shows <10% reducible, report where the cost is irreducible and why, then stop.

## Shared constraints
- Write only inside your assigned directory. Do not read or search the other concurrent directories (`inv-width`, `dqi-dispatch`, `dqi-core`, `chem-audit`). Everything else may be read; `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Start from the round-4 Codex-selected artifacts. Local reference **9.3677221531**: inverse 182,956.5272, DQI 1,598,204,229,385, chemistry 708,547,388.9553. Preserve them until a replacement passes.
- The public leaderboard top is **18.9**, an 8.2× larger product of task ratios than ours, so better solutions exist. A failed approach is evidence about that design only; report it that way.
- Score with the actual objective: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Develop on sampled inputs and short runs, but accept only on full default native verification: DQI needs all 102,400 shots plus independent salts; chemistry needs the native ANGLE/oracle and physical checks. No hash or salt grinding.
- Four agents share 8 cores and one GPU. Respect your worker cap, and wrap full-default verifier runs in `flock /tmp/r6-fullverify.lock …` (the only permitted write outside your directory).
- Report a numeric first checkpoint within a few hours. Maintain `RESULTS.md`: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.