# Round 5: radical paths to 17+

## Bottom line

17.7 looks harder after the round-4 breakdown: under the host cost model chemistry tops out nearer 1.6× than 2×, so 17+ needs DQI near 3× on top of strong gains elsewhere. My range for round 5 is about **12 (conservative) to 14 (target)**, with \~16.7 as a stretch.

- **The most actionable find is public.** Gu, Kemeny, Yelin and Lukin released an Apache-2.0 DQI decoder for this exact Reed–Solomon instance with a division-free root sweep at **23 Toffolis per point** (ours: 88), a one-repair-per-run zero-pivot fix, and measured per-round coefficient windows. At 1,837 qubits it scores worse than ours under Q², so the job is transplanting its ideas into our \~961-qubit layout.
- **DQI correctness comes first.** Cutting 380 failures to ≤8 per 102,400 shots moves R from 1.697 to 1.223 (1.39×). It also removes hash-to-hash noise of about ±1.5% in cost that currently hides smaller gate wins.
- **Chemistry's bottleneck is inner state preparation, not ANGLE.** It is 57% of Toffolis (17,333) and its QROAM temporaries are 337 of 478 qubits. Its cost scales with the W-table size R·C·(B+1), so the lever is a smaller table: B ≤ 511 in one block, or two blocks of \~351 vectors and \~175 columns (about 1.4–1.45× at fixed λ\_eff).
- **Quantization is the binding constraint on every smaller shape.** The inherited basin reaches \~11.3 mHa in float but \~19.6–20 mHa after 8-bit rounding, so \~8–9 mHa of the 20 mHa budget goes to quantization. Cutting that loss is the chemistry enabler; my earlier lookup-pruning bet only touches the 35% ANGLE share and is demoted.
- **Inverse square root:** Horner is 58% of Toffolis in the profiled design and width remains the main lever, worth about 1.1–1.3×.

Suggested round-5 order: DQI exception repair; chemistry quantization-loss reduction on the incumbent; a leaner compact DQI decoder; smaller chemistry shapes once quantization improves; then the width-neutral division-free sweep.

## Score arithmetic

Reaching 17.71 from 9.37 needs the product of the three improvement ratios to grow **6.76×**; 17.0 needs 5.98× and 14.0 needs 3.34×. Because the score is a geometric mean, a 2× gain on any one task is worth exactly as much as a 2× gain on any other.

| Task | Local ratio now | Cost now | Dominant cost factors now |
| --- | --: | --: | --- |
| Inverse square root | 11.35 | 182,957 | Q = 90, \~1,782 Toffolis + \~251 Toffoli-eq of rotations |
| DQI | 15.29 | 1.598 × 10¹² | R = 1.697, Q = 961, T ≈ 1.020 M |
| Chemistry | 4.74 | 708.5 M | λ\_eff = 48.59, T = 30,504, Q = 478 |
| Combined | **9.37** |  | product of ratios = 822 |

One combination that reaches 17.7: DQI 3.25× (ratio \~50), chemistry 1.6× (\~7.6) and inverse square root 1.3× (\~14.8). The round-4 breakdown suggests chemistry tops out near 1.6× under the host cost model, which pushes the burden onto DQI. Each of those is a step change, not a polish; no single task can plausibly supply 6.8× alone.

Hard floors that the formulas themselves impose:

- **DQI repetition.** At default shots, zero observed failures still gives R = 1.19126. The current R of 1.697 leaves **1.42×** on the table from correctness alone.
- **DQI width.** The 840 syndrome qubits are a floor, so Q² cannot fall more than 1.31× from 961. A design that stores Λ and Ω (71 coefficients × 12 bits = 852 qubits) plus a kernel lands near 950, so 961 is already close to its design floor.
- **DQI combined floor.** At Q = 950 and R = 1.19, cost is 1.075 × 10⁶ × T. T = 600 k gives ratio \~38; T = 400 k gives \~57.
- **Chemistry normalization.** Codex's bound of 46.79 sits 3.7% below the current 48.59. Nearly all chemistry headroom is in T × Q, now 14.6 M versus 68.1 M at the host baseline.
- **Inverse square root synthesis.** Rotations are \~12% of cost, so dropping all of them would gain at most \~1.14×. The arithmetic and the 90-qubit width carry the rest.

The practical reading: DQI has the clearest path to a large factor (correctness plus a cheaper phase scan), chemistry needs a smaller inner-preparation table, which first needs lower quantization loss, and inverse square root needs width more than gates.

## Research notes

Three public sources change the picture: Gu et al.'s circuits, the DFTHC cost model in Low et al., and Google's own statement that its internal circuits beat everything published.

**Google's own runs.** The [announcement post](https://opensource.googleblog.com/2026/10/announcing-verifiable-quantum-circuit-compilation-competition-on-kaggle.html) says internal Gemini runs improved all three benchmarks beyond published circuits, with a hash of those circuits in the [companion paper](https://arxiv.org/abs/2610.23832) to be revealed later. Treat 17.7 as plausible for a strong team, not as evidence of an exploit. The paper also confirms classical bit operations and a condition stack are part of the kickmix gate set, and that the private salt exists specifically to block seed grinding.

**Gu et al. (Harvard and Anthropic).** [arXiv:2610.22367](https://arxiv.org/abs/2610.22367) and the [andigu/uncompute](https://github.com/andigu/uncompute) repo cover eight tasks. Their three recurring design principles map directly onto our bottlenecks:

1. Conditional routing: move branch differences into cheap controlled swaps so one uncontrolled arithmetic block serves every branch.
2. Decode branch history from the evolving state during uncomputation instead of storing it.
3. Spend the failure budget by first localizing rare-case resources (extra iterations, stack entries, lanes), then omitting them.

**Their DQI decoder** ([task](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/README.md), [circuit notes](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md)) is the same RS(4095, 4025) instance and weight law. Stage costs, Toffolis: Berlekamp–Massey 141,071 forward and 86,176 reverse; in-place ω transform 15,246 each way; half-degree Euclid for Ψ 44,470 each way; sweep 92,552 (4,024 points × 23). Peak is 1,837 qubits, set during the Euclid passes. Field multiply is 42 Toffolis and inversion 66. Failure sources are measured coefficient cutoffs (3.8e-5), a one-entry discrepancy stack (1.0e-4), and a second zero pivot (4.5e-6). Under their grader the weight law has support 19–35, mode 32 at 18.4%, weight 35 at 2.8%.

Two differences matter when porting. Their mask is a fixed constant, so they skip 71 points and fold constants at build time; Kaggle masks are per-shot classical registers. And their width is \~2× ours, because σ, the BM template and Ψ's Bézout pair all live alongside the 840 syndrome qubits.

**Their other tasks.** The [Coulomb task](https://github.com/andigu/uncompute/blob/main/tasks/coulomb-potential/README.md) (1,140 Toffolis, 142 qubits) uses 128 lines with a 6-bit local coordinate, but its relative tolerance is 32–14,000× looser than ours, so its table does not transfer. The [Fe₂S₂ chemistry task](https://github.com/andigu/uncompute/blob/main/tasks/chemistry-fe2s2/README.md) cut an alias lookup from 41,050 to 10,271 Toffolis by scheduling rows so aligned address subtrees repeat, and by keying cheap Clifford "frames" off address bits. Both tricks are relevant to our ANGLE oracles.

**Low et al. DFTHC cost model** ([arXiv:2502.15882](https://arxiv.org/abs/2502.15882)). Rotation lookup costs about N + RB per call; inner preparation uses QROAM costing RC⌈(B+1)/2^k₂⌉ + 2^k₂·b₂ with (2^k₂ − 1)·b₂ temporary qubits; Givens rotations cost 2(N−1)·b\_rot per call. The round-4 report confirms the host follows this: inner preparation is 57% of our Toffolis and its QROAM temporaries are 337 of our 478 scored qubits. Their FeMoCo-76 solution, (R, B, C) = (15, 57, 19), reports Λ = 179.7, E\_gap = 5.38 and λ\_eff = 43.65. Those numbers are not directly comparable (optimized symmetry shifts, a different accuracy metric), but our implied Λ of \~303.5 is far above 179.7, so the assumptions behind Codex's 46.79 bound deserve a check.

## DQI: radical ideas

DQI can plausibly reach a 1.9–2.9× task gain (ratio 29–44) at ≤\~960 qubits by stacking four changes. Iteration 4 already sits near the width floor for a design that stores σ and ω (852 qubits), so every idea below must be width-neutral.

Width sensitivity at Q ≈ 961: +10 qubits costs 2.1%, +100 costs 21.9%. Repetition factor versus failures in 102,400 shots (values for 20, 50 and 100 are my Poisson approximation):

| Failures | 0 | 8 | 20 | 50 | 100 | 380 |
| --- | --: | --: | --: | --: | --: | --: |
| R | 1.191 | 1.223 | 1.255 | 1.31 | 1.39 | 1.697 |

### D1. Exception repair inside the compact decoder

The 380 phase errors come from the degree-one-quotient assumption. Gu et al. handle a zero pivot (0.36% of shots in their half-degree Euclid) with one controlled realignment per run plus a six-bit record of where it happened, leaving 4.5e-6 for a second pivot. Port that pattern: detect the first leading-coefficient defect, realign once, record the position, replay it in reverse.

- Target: ≤8 failures, ≤3% more Toffolis, ≤10 more qubits.
- Gain: R 1.697 → 1.223 gives 1.34–1.39× on DQI, about 9.37 → 10.4 overall.
- Kill: more than \~40 extra qubits or 5% extra Toffolis.

D1 also removes measurement noise. With \~380 expected failures the count varies by about ±20 between circuit hashes, which moves R by about ±0.026 and cost by about ±1.5%. That is why Codex's 87-Toffoli kernel scored worse than the 88-Toffoli one; below \~10 failures R is stable at 1.19–1.23 and small gate wins become visible.

### D2. A leaner compact decoder

The decoder (≈1.02M minus the \~360k sweep, so \~660k forward plus reverse) is the largest block. Gu's BM costs 227k forward plus reverse, and our own iteration-3 BM spent 293.5k (159.7k forward, 133.8k reverse, per Claude's profile), but both at much greater width. Port the width-independent tricks into the 961-qubit decoder:

- **Measured coefficient windows.** Simulate the classical decoder under the weight law and clip each step's active coefficient range where the clipped mass is below \~4e-5. Weights below 19 essentially never occur, so early steps are generic.
- **Shared arithmetic path.** Route branch differences through controlled swaps so one inversion, one multiply and one multiply-add serve both branches.
- **History from state.** Gu recover a displaced discrepancy from the invariant σ₀·δ = 1 instead of storing it; look for the analogous invariant in our Euclid layout.
- **Measurement-erase forward temporaries** and recompute only those in the reverse pass.

Target: decoder ≤420k at unchanged width, roughly 1.3× on DQI.

### D3. The sweep

**(a) Per-point cost 88 → 80.** Iteration 4 pays 88 per point; Gu's full per-point cost with a division is 80 (11 for the root test plus 69 for division and trace pairing). Closing that gap saves \~33k. An AND-count search (SAT or ABC) on the per-point function Tr(ν·ω/σ') might reach \~55–60.

Claude's profile splits the inherited 80 into an 11-Toffoli zero test, 66 for the inverse-parity computation and 3 for the final parity, the same split as Gu's division path. The compact layout started at 92 and Codex's GF(16)-cubic trace kernel brought it to 88; first establish whether the compact storage basis forces that extra cost.

**(b) Division-free sweep, width-neutral.** Write σ = E² + zO², so σ' = O² and the halves have degree ≤17. A Euclid pass on (E, O) gives Ψ with Ψ·σ' = 1 at every root, costing \~44k per direction in Gu's circuit. The sweep then needs no division: 11 for the root test plus 12 gated pairings, 23 per point (\~94k instead of \~330–360k).

The catch is width: Gu hold Ψ's Bézout pair (36 words, 432 qubits) in syndrome words freed by their in-place σS transform, and our compact layout has no such free space. The radical variant: **replace ω by F = ω·Ψ mod σ in place**. Multiplying by Ψ is invertible in GF\[z\]/σ, so F fits in ω's 35 words. The sweep then needs only σ and F, and the phase becomes a classically controlled linear parity of F(x) gated by the root flag: about **11 Toffolis per point**. The open problem is doing that transform with fewer than \~100 extra peak qubits (break-even: peak width may grow by about √(T\_old/T\_new)).

Round 4 reached the same width conclusion from two directions: Claude found the upstream sweep needs extra symbol storage that the Q² penalty erases, and AGY's break-even for it was \~1,723 qubits against the 1,429-qubit design. Only a variant that adds almost no peak width is worth building.

### D4. Width trimming

The scan needs \~880 qubits (σ, ω and a small kernel) but the circuit peaks at 961. A liveness profile should show which phase sets the peak; each 10 qubits removed is worth 2.1%.

### Projection

| Scenario | Q | T | R | Cost (×10¹²) | DQI ratio | Task gain |
| --- | --: | --: | --: | --: | --: | --: |
| Now | 961 | 1.02M | 1.697 | 1.598 | 15.3 | 1.00× |
| D1 only | 965 | 1.05M | 1.223 | 1.196 | 20.4 | 1.34× |
| D1 + D2 + 80 per point | 950 | 0.77M | 1.223 | 0.850 | 28.8 | 1.88× |
| + width-neutral F sweep | 950 | 0.61M | 1.223 | 0.673 | 36.3 | 2.37× |
| Stretch | 930 | 0.53M | 1.21 | 0.555 | 44.0 | 2.88× |

Porting notes: Kaggle masks are per-shot classical registers, so Gu's 71-point skip and build-time constant folding become classically conditioned gates. Confirm the field polynomial and position indexing match the host verifier before reusing any Gu component. Keep Apache-2.0 attribution, as iteration 3 did.

## Chemistry: radical ideas

Inner state preparation dominates chemistry cost, its size is fixed by the factorization shape, and every cheaper shape tried so far was blocked by \~8–9 mHa of 8-bit quantization loss. The plan is therefore: cut quantization loss first, then shrink the W table. Normalization is not a lever (gap 3.915157 against a 3.914898 floor; λ\_eff within 4% of the fixed-shift bound).

### Where the cost is

Selected candidate R1/B719/C350, 8-bit angles, J = 7, from the native verifier output:

| Component | Toffolis | Share | Persistent qubits | Temporary qubits |
| --- | --: | --: | --: | --: |
| Inner state preparation | 17,333 | 57% | 45 | **337** (sets the maximum) |
| ANGLE, forward + reverse | 10,604 | 35% | 56 | 57 |
| Givens rotations | 2,250 | 7% | — | — |
| Outer state preparation | 317 | 1% | 40 | 90 |
| Total | 30,504 |  | 141 | 337 → Q = 478 |

How the host prices inner preparation, reconstructed from the report's formulas (it reproduces 17,333 exactly):

- Cost ≈ R·C·⌈(B+1)/2^k₂⌉ + 2^k₂·b₂ + cleanup (R·C + B + 1 + b₂ + N) + 8·⌈log₂(B+1)⌉. For the incumbent the first term is 350 × 45 = 15,750, over half of all Toffolis.
- The host takes the largest k₂ ≤ 7 with (2^k₂ − 1)·b₂ ≤ b\_rot·(N−1) = 600, where b₂ = ⌈log₂B⌉ + 11. With b₂ = 20–21 that gives k₂ = 4 and \~300–315 temporary qubits. There is no submission knob.
- Two width cliffs follow. Nine-bit angles raise the threshold to 675 and allow k₂ = 5 (Q ≈ 820). Fewer than 257 vectors per ensemble give b₂ = 19, which also allows k₂ = 5. Stay at ≤8 bits and ≥257 vectors per ensemble.
- So width is pinned near 460–480 for this architecture, and gains must come from T, mainly the R·C·(B+1) table.

### What smaller shapes would buy

My estimates from those formulas, assuming ANGLE Toffolis scale with address count and λ\_eff is unchanged; accuracy is untested:

| Shape | Inner prep | ANGLE | T | Q | T×Q vs now |
| --- | --: | --: | --: | --: | --: |
| R1/B719/C350, 8-bit (now) | 17,333 | 10,604 | 30,504 | 478 | 1.00× |
| R1/B511/C350, 8-bit | \~12,550 | \~7,830 | \~22,950 | \~461 | 0.73× |
| R2/B351/C175, 8-bit | \~8,890 | \~10,380 | \~21,840 | \~461 | 0.69× |
| R2/B351/C175, 7-bit | \~8,890 | \~9,080 | \~20,240 | \~461 | 0.64× |
| R4/B180/C88, 8-bit | k₂ = 5 cliff |  |  | \~740 | worse |

These match Claude's earlier 1.38–1.49× predictions for the same shapes. Column-blocked R2 shapes keep every vector free to span all 76 orbitals; only each column is confined to one block.

### C1. Cut quantization loss (the enabler)

Relaxing the inherited solution in float reaches 11.3 mHa, but the best 8-bit version sits at \~19.6–20.0 mHa, and Claude's quantizer gave 22.7 from the same basin. Two-block float fits at 12.9 and 14.8 mHa quantized to 25 and 28. Target: on the incumbent shape, quantized error within \~4 mHa of float, so float fits up to \~16 mHa become submittable.

1. **Error-feedback rounding along the Givens tree.** Round angles from the root outward and re-solve each vector's remaining continuous angles after every rounding, so later angles absorb earlier errors.
2. **Progressive freezing.** Quantize a fraction of vectors, refit weights, offsets and the still-continuous vectors, and repeat.
3. **Optimize the tree itself.** Orbital order and tree topology are free and shared by all addresses; choose them for low quantization sensitivity.
4. Score every step on the weighted two-body residual, never on per-angle or per-vector error.

Track "quantization loss" (quantized minus float error) as a first-class metric. Kill if the incumbent's loss stays above \~6 mHa.

### C2. Shrink the W table once C1 works

- Fit R1/B511 and R2/B351/C175, plus neighbours on the safe side of both cliffs, from warm starts. Prune slowly with leverage scoring (Claude's best: 12.9 mHa float at B479), then apply the improved quantizer.
- Build the full host-model frontier over (R, B, C, b\_rot, J) with the current signed-gauge oracle costs and the host's own QROAM choices; the report notes none exists yet.
- Watch λ\_eff: removing only 3–5 columns raised it 1.5–2.3% in AGY's C347/C345 runs.
- Expected 1.35–1.55× on the task if λ\_eff holds.

### C3. Seven-bit angles after C1

Seven bits cut Givens to \~1,950 Toffolis and ANGLE by about one eighth (\~7% of T) without moving k₂. Worth it only once quantization loss is under control.

### Demoted

- **Locality-structured dictionary with pruned lookups.** It targets only the 35% ANGLE share, which already runs at \~482 Toffolis per oracle for 795 addresses, and it conflicts with block structure if blocks follow orbital clusters.
- **Select-swap inside the 128-qubit cap and address reordering.** The signed-gauge oracles already capture most of this.
- **Normalization.** Low et al.'s Λ of 179.7 relies on optimized symmetry shifts that the verifier fixes.

Chemistry projection: conservative 1.15×, target \~1.35×, stretch \~1.6× (R2/B351/C175 at 7 bits with λ\_eff unchanged).

## Inverse square root: radical ideas

Expect 1.1–1.3× here. At 90 qubits and \~1,782 Toffolis the arithmetic is close to saturated, so width is the main lever: removing 10 qubits (11%) pays for about 220 extra Toffoli-equivalents.

Profiled stage costs of the inherited 94-qubit design: normalization 210, coefficient lookup 308, Horner 1,043 (58%), output shift 228, plus 34 rotations. Horner is the gate target. AGY's narrower 86–90-qubit candidates cost \~2,215 Toffoli-equivalents (\~199k overall, worse than the incumbent), while Codex's 90-qubit circuit won through a better storage schedule; width must be bought with scheduling, not with extra gates.

### I1. Width-first redesign

Profile which registers are live at the 90-qubit peak, then redesign that phase rather than the whole circuit. Gu et al. ship liveness and occupancy profilers that are worth copying.

- Normalize x in place and restore it, so the offset t is a slice of x rather than a copy.
- Accumulate Horner terms into the looked-up coefficient registers in place.
- Where carry ancillas set the peak, test low-ancilla (Cuccaro-style) adders: about n more Toffolis per addition in exchange for about n fewer qubits.

### I2. Power-of-two-aligned, non-uniform intervals

The tolerance is relative, so the allowed interval width grows roughly in proportion to the mantissa. Uniform intervals within an octave therefore over-provision its upper end. Refining intervals as a binary prefix tree (as Gu's HHL-angle circuit does) could cut table rows by 20–30% and narrow coefficient ranges; lookup cost is about one Toffoli per row.

### I3. Exploit mean-executed accounting

Fixups conditioned on random measurement outcomes run half the time on average. Prefer measurement-based erasure plus fixup over explicit reversal wherever the fixup costs less than twice the reversal.

### What not to chase

- Rotation synthesis is \~251 Toffoli-equivalents (12%). Doubling the synthesis error budget saves only \~0.14 per rotation, about 5 in total.
- Gu's Coulomb circuit (128 lines, 1,140 Toffolis) meets a much looser relative tolerance over an 18-bit domain; its line table does not meet 2^-24 relative accuracy over 22 bits.
- Goldschmidt, Walsh and Boolean-synthesis families were already tried without a competitive complete circuit.

## Ideas from rounds 3–4 that were not fairly evaluated

Most structured chemistry fits failed at quantization rather than in float; one (orbital-disjoint blocks) failed for a genuine structural reason and should not be retried.

| Idea | Reported result | What it does and does not show | Fair re-test |
| --- | --- | --- | --- |
| AGY disjoint-block chemistry fit | Stalled at 479.36 mHa | It is evidence, but only against that partition: disjoint blocks drop the cross-block correlations, so inter-block Coulomb terms cannot be represented. Column-blocked R2 shapes keep vectors spanning all orbitals | Do not retry; test column-blocked R2/B351/C175 instead (chemistry C2) |
| Claude two-block fits | 13.0 mHa float, 25.0/27.8 mHa after quantization | The same local quantizer also degraded an inherited control from \~19.6 to 22.7 mHa | Run through the round-3 integer-angle pipeline |
| Grouped dictionary, 28 groups | 17.02 mHa float, 391 mHa at 7-bit rounding | Direct rounding at 7 bits, and covariance refits also missed the cap. R24/B64 also doubles ANGLE addresses to 1,612, so it may lose even if accurate | Low priority; prefer \~700-address R2 shapes after quantization improves |
| AGY rational-reconstruction DQI | Projected 938–1,312 qubits and scores of 14–17 | Unbuilt; the overall-score arithmetic was also off | Compare against the D1–D3 path before building |
| Codex Wronskian padding (DQI) | Adapter 28,069 → 7,140 Toffolis per direction | Component only; no complete decoder | Finish only if it enables the width-neutral F sweep |

A fair test protocol for any structured chemistry idea:

1. Start from the current solution projected onto the structure, not from scratch.
2. Quantize with integer-angle moves plus weight and offset re-solves, as in round 3.
3. Score only with the native verifier on quantized words.
4. Compare at equal λ\_eff × T × Q from the host model, not at equal shape.
5. Report float, rounded and refit errors side by side, so a quantizer failure is visible as one.

## Round-5 Codex briefs

Seven targeted briefs, re-ranked after the round-4 breakdown. Overall-score effects assume the other two tasks stay fixed at 9.37.

| # | Task | Goal | Task gain | Overall | First checkpoint | Kill criterion |
| --- | --- | --- | --: | --: | --- | --- |
| B1 | DQI | Exception repair in the compact decoder: ≤8 failures at Q ≤ 970, T ≤ 1.05M | 1.34× | →10.3 | Failure census by cause on ≥1M law-sampled shots | >40 extra qubits or >5% extra T |
| B2 | Chem | Cut 8-bit quantization loss on the incumbent to ≤4 mHa over float (error feedback, progressive freezing, tree choice) | enabler | — | Loss of the current quantizer vs each new method on the same float basin | Loss stays >6 mHa |
| B3 | DQI | Leaner compact decoder: forward + reverse ≤420k at unchanged width | \~1.3× | +9% | Per-step Toffoli and live-qubit table; measured coefficient windows | <10% decoder saving after windows + routing |
| B4 | Chem | Smaller W table: R1/B511 and R2/B351/C175 with the B2 quantizer; full host-model frontier with current oracle costs | 1.35–1.55× | +10–16% | Float fits ≤16 mHa at frontier shapes | Quantized >20 mHa at every frontier shape |
| B5 | DQI | Per-point cost ≤80; research the width-neutral F = ω·Ψ mod σ sweep | 1.03–1.3× | +1–9% | Peak-width plan for the F transform on paper | Transform needs >\~110 extra peak qubits |
| B6 | Inverse | Width-first schedule plus cheaper Horner (58% of Toffolis) and power-of-two-aligned intervals | 1.1–1.3× | +3–9% | Liveness profile at the 90-qubit peak | Cannot reach Q ≤ 82 at T ≤ 2,000 |
| B7 | Chem | Seven-bit angles once B2 works | \~1.07× | +2% | Quantization loss at 7 bits on the B4 winner | Loss pushes error over 20 mHa |

Scenarios for the combined score:

| Scenario | DQI | Chemistry | Inverse | Score |
| --- | --: | --: | --: | --: |
| Conservative (B1, partial B3, partial B6) | 1.6× | 1.15× | 1.1× | 11.9 |
| Target (B1, B3, per-point 80, B2 + one B4 shape) | 2.1× | 1.35× | 1.2× | 14.1 |
| Stretch (width-neutral F sweep, R2 at 7 bits) | 2.7× | 1.6× | 1.3× | 16.7 |

Reaching 17.7 from the stretch row needs DQI at about 3.25× instead of 2.7×: roughly 470k Toffolis at \~930 qubits with ≤8 failures, or a chemistry gain outside the host model's shape space.

Prompting notes for each brief:

- Give the full objective and the break-even rule (for DQI, a width increase must be paid for by T falling by at least (Q\_new/Q\_old)²).
- Require a first checkpoint within a few hours, reported as numbers, before any long build.
- For chemistry, track float error, quantized error and their difference for every candidate.
- Point at the relevant Gu et al. README sections and require attribution for anything reused.
- Report component wins and end-to-end wins separately, so a component result is never quoted as a score.
- Run B1/B3 in one DQI workspace, B5 in a second, chemistry (B2 → B4 → B7) in a third, and B6 in a fourth.

## Evaluation guardrails

Round 4's dead ends mostly traced to proxies (float errors, component counts, projections) standing in for the full objective; these rules keep round 5 honest.

- **Score with the full objective every time:** Q × (T + T-count/4), R × Q² × T, and λ\_eff × T × Q. Report a break-even for any width increase.
- **Estimate DQI failures before the final run.** Use ≥1M law-sampled local shots so the expected failure count, and hence R, is known before the 102,400-shot verification. Never infer rates from small instances or short runs.
- **Width claims need a liveness profile,** because the scored Q is the circuit's total qubit count, set by its peak.
- **Chemistry evidence is native quantized verification only.** Run the integer-aware refit before declaring any structure dead.
- **No seed grinding.** Official salts are private, so local luck does not transfer.
- **Keep anchor tests:** the zero syndrome, single-symbol errors, and the weight-35 tail for DQI; extreme x and interval boundaries for the inverse square root.
- **Keep attribution** for anything adapted from Apache-2.0 sources.

## Sources

- Khattar et al., [Leveraging verifiable structure in quantum circuits to find more efficient compilations with agentic AI](https://arxiv.org/abs/2610.23832), arXiv:2610.23832
- Gu, Kemeny, Yelin and Lukin, [Discovering state-of-the-art quantum circuits with large language models](https://arxiv.org/abs/2610.22367), arXiv:2610.22367
- [andigu/uncompute repository](https://github.com/andigu/uncompute)
- [dqi-reed-solomon task README](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/README.md) and [discovered-circuit notes](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md)
- [coulomb-potential task README](https://github.com/andigu/uncompute/blob/main/tasks/coulomb-potential/README.md)
- [chemistry-fe2s2 task README](https://github.com/andigu/uncompute/blob/main/tasks/chemistry-fe2s2/README.md)
- Low et al., [Fast quantum simulation of electronic structure by spectrum amplification](https://arxiv.org/abs/2502.15882), arXiv:2502.15882
- Khattar et al., [Verifiable quantum advantage via optimized DQI circuits](https://arxiv.org/abs/2510.10967), arXiv:2510.10967
- Google Open Source Blog, [Announcing the verifiable quantum circuit compilation competition on Kaggle](https://opensource.googleblog.com/2026/10/announcing-verifiable-quantum-circuit-compilation-competition-on-kaggle.html)
- Round 1–4 reports, host-baseline, profiling and verifier notes supplied in chat