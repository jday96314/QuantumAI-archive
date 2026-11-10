# Round 5 / inv-width — width-first inverse square root

You own `initial_experiments/5/inv-width`.

**Starting facts:** the round-4 Codex-selected inverse circuit uses Q=90, ~1,782 mean Toffolis and 34 rotations; cost 182,956.5272 (max error/tolerance ≈0.9797). Cost = Q × (mean_Toffolis + total_T/4), so cutting 10 qubits (11%) pays for ~220 extra Toffoli-equivalents. An inherited 94-qubit design profiled at normalization 210, lookup 308, Horner 1,043, output shift 228 Toffolis. No agent has attacked this task's width directly since round 4.

## Work plan
1. Reproduce the selected cost. Produce a liveness profile showing exactly which registers are live at the 90-qubit peak, plus per-stage Toffolis.
2. Width probes on the incumbent:
   - **Segmented phasing:** carries flow only upward, so phase a low segment of the result in a narrow register, keep its carry-out, uncompute it, then compute and phase the high segment with that carry in.
   - Low-ancilla (Cuccaro-style) adders wherever carry ancillas set the peak.
   - In-place normalization and Horner.
   - Measurement-based erasure with outcome-conditioned fixups (only executed Toffolis count).
3. Seed at least three algorithm families **without the incumbent's code**, for example: direct quadratic with segmented phasing; bit-serial evaluation on narrow registers; table-heavy with power-of-two-aligned non-uniform intervals and measurement-based erasure. Fit coefficients against each design's exact integer pipeline. Score every family on the full objective at its first checkpoint, and cross-pollinate ideas between them.
4. Develop with exhaustive integer-model checks and sampled-domain native runs. Run the full verification (200,540,928 shots) only for designs projected to beat the incumbent.

## Success and scope
- Target Q ≤ 70 at ≤2,000 Toffoli-equivalents (≈1.3×). Stretch Q ≤ 60 (≈1.5×).
- If, after the liveness pass and three families, nothing projects below ~165k, document a concrete obstruction for each family and stop.

## Shared constraints
- Write only inside your assigned directory. Do not read or search the other concurrent directories (`inv-width`, `dqi-dispatch`, `dqi-core`, `chem-audit`). Everything else may be read; `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Start from the round-4 Codex-selected artifacts. Local reference **9.3677221531**: inverse 182,956.5272, DQI 1,598,204,229,385, chemistry 708,547,388.9553. Preserve them until a replacement passes.
- The public leaderboard top is **18.9**, an 8.2× larger product of task ratios than ours, so better solutions exist. A failed approach is evidence about that design only; report it that way.
- Score with the actual objective: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Develop on sampled inputs and short runs, but accept only on full default native verification: DQI needs all 102,400 shots plus independent salts; chemistry needs the native ANGLE/oracle and physical checks. No hash or salt grinding.
- Four agents share 8 cores and one GPU. Respect your worker cap, and wrap full-default verifier runs in `flock /tmp/r6-fullverify.lock …` (the only permitted write outside your directory).
- Report a numeric first checkpoint within a few hours. Maintain `RESULTS.md`: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.