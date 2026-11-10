# Round 5 / dqi-dispatch — exception-correct compact DQI decoder

You own `initial_experiments/5/dqi-dispatch` and are the only agent modifying a copy of the compact decoder.

**Starting facts:** the selected decoder is `initial_experiments/4/codex/build_cubic_compact_decoder.py` (core `compact_product_decoder.py`). It has Q=961, mean T≈1,019,866, 380 phase / 0 value errors, R=1.69685, cost 1.598204229385e12.

`initial_experiments/5/dqi-repair` diagnosed the failures (read its RESULTS.md and code):
- The 380 phase errors come from 744 wrong decodes; 99.7% are an isolated degree-two quotient.
- Its degree-(2,1) block (`gap_cofactor.py`, `gap_block.py`, `classical_repair.py`, `census.cpp`) cuts failures to 2 in a classical splice.
- Its forced stage-10 circuit fixed targeted exceptions but broke ordinary inputs, and was 10,506 T cheaper than the regular block it replaced.
- Unresolved: dispatch at any stage, immediate completion after a degree-two quotient, the final two slots, block-straddling gaps, and the serial region (k ≥ 22).

Width facts from `initial_experiments/5/dqi-hybrid`:
- Peak 961 = 840 + 22 completion flags + 99 temporaries (including a 72-bit block head) at block k=19.
- The serial region stays ≤914.
- The scan holds 959 live qubits: 840 + 35 flags + 84 scratch.

## Work plan
1. **Uniform local-transfer block.** Every block computes its whole 2×2 transfer from local moments, using a small reversible local Euclid that routes the (1,1,1), (2,1), (1,2) and (3) quotient patterns through one arithmetic path, then applies it with the shared cofactor passes. No branch between regular and gap blocks.
   - Handle block-straddling gaps with a carry bit plus a controlled one-word realignment.
   - Use two-slot look-ahead steps in the serial region.
   - If quantum-selected transfer between the regular and (2,1) blocks proves cheaper, use that instead.
2. Build a classical model of the design first. It must reproduce ~2 failures on the exact 102,400 host shots and ≤~20 per million law-sampled cases.
3. **Width trim**, co-designed with step 1: recompute completion from the polynomial state (a zero test per step) or keep a 6-bit counter instead of 35 unary flags, and shrink the block head near the peak.
4. Emit native components, then the full circuit. Check zero-syndrome, one-error, weight-35 and ordinary controls with independent random masks, then run the full native verification.

**Budget:** at Q=961 and ~2 failures, ≈+47k T still gives 1.35× and ≈+88k T gives 1.3×. Each 10 qubits removed is worth ~2.1%.

## Success and scope
- Target ≤10 failures, Q ≤ 965, T ≤ 1.09M (≥1.3×). Stretch ≥1.45× with the width trim.
- Kill if the classical model cannot get below ~20 failures on the host shots, or the projected overhead exceeds 90k T at unchanged width.
- The user may relay an attribution table from a parallel read-only core analysis. Do not read other new directories.

## Shared constraints
- Write only inside your assigned directory. Do not read or search the other concurrent directories (`inv-width`, `dqi-dispatch`, `dqi-core`, `chem-audit`). Everything else may be read; `from_host/` and prior artifacts are read-only; never modify the verifier, confidence or scoring.
- Start from the round-4 Codex-selected artifacts. Local reference **9.3677221531**: inverse 182,956.5272, DQI 1,598,204,229,385, chemistry 708,547,388.9553. Preserve them until a replacement passes.
- The public leaderboard top is **18.9**, an 8.2× larger product of task ratios than ours, so better solutions exist. A failed approach is evidence about that design only; report it that way.
- Score with the actual objective: DQI `R(epsilon_plus) * Q**2 * mean_Toffolis` (host Clopper–Pearson); chemistry `lambda_eff * T_total * Q_total`; inverse `Q * (mean_Toffolis + total_T/4)`.
- Develop on sampled inputs and short runs, but accept only on full default native verification: DQI needs all 102,400 shots plus independent salts; chemistry needs the native ANGLE/oracle and physical checks. No hash or salt grinding.
- Four agents share 8 cores and one GPU. Respect your worker cap, and wrap full-default verifier runs in `flock /tmp/r6-fullverify.lock …` (the only permitted write outside your directory).
- Report a numeric first checkpoint within a few hours. Maintain `RESULTS.md`: hypotheses, changes, checks, measured versus estimated numbers, failures and reproduction commands. Attribute reused Apache-2.0 code.