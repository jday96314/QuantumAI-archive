# Round 4 optimization report

## Outcome and scope

The three agent workspaces (`agy`, `claude`, and `codex`) started from the round-3 submission scoring **9.0096229084**. Only Codex produced a better verified submission: **9.3677221531**, a **3.97%** overall increase. The original 14.0 target was not reached; Claude later cleared that goal and pursued 9.4, which was also not reached. AGY and Claude preserved the inherited submission rather than replacing it with a worse or invalid candidate.

The score is the geometric mean of the three reference-to-submission cost ratios. Lower task costs are better; gate reductions alone do not establish a score improvement. Inverse square root charges qubits times Toffoli-equivalent work, DQI charges repetition factor times qubits squared times Toffolis, and chemistry charges effective normalization times Toffolis times qubits.

| Agent | Best packaged score | Inverse-square-root cost | DQI cost | Chemistry cost |
| --- | ---: | ---: | ---: | ---: |
| Inherited / Claude / AGY | 9.0096229084 | 190,711.9555 | 1.612737419463e12 | 757,163,489.50 |
| Codex | **9.3677221531** | **182,956.5272** | **1.598204229385e12** | **708,547,388.96** |

These results come from the local [Codex summary](codex/submission.summary.json), [Claude summary](claude/submission.summary.json), and [AGY score summary](agy/logs/submission_score_summary.json). This report synthesizes recorded experiments; it does not rerun the native verification campaigns.

## Codex: improvements across all three tasks

### What worked

- **Inverse square root: affine lookup compression and workspace reuse.** Subtracting a cheap affine function of the interval address reduced the constant coefficient from 27 to 23 bits. Delayed coefficient loading and early curvature-bit release reduced width from 94 to 90 qubits. An exhaustive integer-model invariant permitted shorter modular offset restoration; measurement-based cleanup reduced output-shift work. The selected circuit uses 34 rotations and about 1,782 mean Toffolis, passes **200,540,928 native shots**, and cuts task cost **4.07%**. Its maximum error/tolerance ratio is about 0.97967.
- **DQI: compact Euclidean decoding and a cubic-extension phase kernel.** Consuming syndrome words into compact cofactors, blocking early Euclidean steps, switching to serial updates late, and deferring completion handling reduced width from 1,429 to **961 qubits**. A GF(4096)-over-GF(16) trace construction reduced the scan kernel from 92 to 88 Toffolis per field position. The full decoder passes **102,400 native shots**, with zero value errors and **380 phase errors** allowed and charged by the scorer. About 1,019,866 mean Toffolis and a repetition factor of 1.69685 consume most of the width savings: task cost improves only **0.90%**.
- **Chemistry: signed angle gauges, lookup batching, and normalization-aware refits.** Equivalent signed angles make two output bits dependent, allowing larger lookup subtrees under the per-oracle 128-qubit limit. Eleven seven-word-batch oracles use 123 qubits each; total scored preparation width falls from 526 to **478**. Integer dictionary updates and factor/offset refits reduce effective normalization to **48.59417**. Despite increasing Toffolis to 30,504, task cost falls **6.42%**. The submission passes **50,880 native shots**; one-body error is 11.9863 mHa and two-body error is 19.9953 mHa, very close to the 20 mHa cap. Each oracle also passes 2,048 random incoming-target checks.

The selected inverse, decoder, and chemistry exporter were independently rebuilt with identical outputs. The final archive and independent audit are documented in [EXPERIMENT_REPORT.md](codex/EXPERIMENT_REPORT.md); the archive SHA-256 is `9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68`.

### Current chemistry resource breakdown and existing shape scans

The authoritative breakdown is already in the selected submission's [full native verifier output](codex/results/norm_cycle2_signed_verify.json), under `resources`. This is the **708,547,388.96-cost selected candidate**, not the earlier signed-gauge checkpoint described at the start of some research notes. Its shape is **R=1, B=719, C=350, brot=8, J=7**, with eleven oracles and fixed keep-value precision `chi=10`.

| Component | Toffolis in the scored full preparation | Persistent qubits | Temporary qubits |
| --- | ---: | ---: | ---: |
| ANGLE family, forward | 5,302 | 56 shared word qubits | 57 |
| ANGLE family, reverse | 5,302 | Same 56; do not add twice | Same 57 |
| Outer state preparation, including cleanup | 317 | 40 | 90 |
| Inner state preparation, including cleanup | 17,333 | 45 | **337** |
| Givens rotations, both passes | 2,250 | Accounted for by the host model's storage convention | No separate additive workspace term |
| **Total** | **30,504** | **141** | **337 shared maximum** |

Thus `T = 10,604 + 317 + 17,333 + 2,250`, and `Q = (56 + 40 + 45) + max(57, 90, 337) = 478`. **Inner preparation sets the temporary-width bottleneck**. Each actual ANGLE oracle uses 123 qubits: 10 selection + 56 word + 57 temporary. The selection register is supplied by preparation and is not added again in the scored total. The 337 temporary qubits are about 70.5% of the total width; reducing ANGLE scratch alone does not reduce that maximum.

The [host preparation model](../../from_host/verifiers/dfthcsos/scoreinfra/costs_prepare.py) makes the QROAM choices; [heuristic_costs and resource_costs](../../from_host/verifiers/dfthcsos/verify.py) connect it to scoring. Substituting the selected shape into that code gives:

- Inner word width `b2 = ceil(log2(B)) + chi + 1 = 21`.
- Unconstrained gate-minimizing `QR(R*C*B, b2)` starts at exponent **k2=7**. The host defaults to `optimize_qubits=True` and lowers this to **k2=4**, a **16-entry block**, to fit its rotation-storage threshold `brot*(N−1) = 600`. It is therefore not a pure Toffoli-minimizing choice.
- Inner temporary storage is exactly `max(0, ceil(log2(B+1))−2) + ceil(log2(R*C*ceil((B+1)/2**k2))) + (2**k2−1)*b2 = 8 + 14 + 315 = 337`. The 315-qubit term is the difference-bank contribution identified in the research question; address/uniform-preparation terms supply the rest.
- Outer forward QROM uses exponent **2** (block size 4); outer inverse uses exponent **4** (block size 16). Inner cleanup's `QR(R*C, B+1)` uses exponent **0** (block size 1). These exponents are derived from the linked source, rather than fields emitted in the native JSON.

Existing resource-only searches and outputs are also available:

| Evidence | What it covers / limitation |
| --- | --- |
| [Claude costscan.py](claude/costscan.py) | Scans R, total dictionary size, total column count, brot, and J; ranks by T×Q without physical accuracy. Uses an approximate unary lookup cost and assumed eight scratch qubits, not the selected signed-gauge oracle's costs. |
| [AGY cost-table script](agy/scripts/dfthc_cost_table.py) and [saved table](agy/logs/dfthc_cost_table.log) | Selected shapes at 10–12 bits, enumerating admissible J under an ordinary-QROM width estimate; reports T, Q, and preparation persistent/temporary splits. Assumes a fixed normalization for illustrative ratios. |
| [AGY multirank scan](agy/scripts/test_multirank_cost_table.py) | Wider R/B/C search at eight bits, restricting total columns to at least 290 and total vectors to at most 800. J is fixed at 13, and its assumed angle layout is not a measured candidate. |
| [Codex grouped-shape resource output](codex/results/grouped_shape_heuristic_costs.json) and [expanded output](codex/results/grouped_expanded_shape_heuristic_costs.json) | Saved host preparation breakdowns over grouped shapes and precision. These are **heuristic-only** costs, with angle-row counts, rather than complete ANGLE-inclusive T/Q predictions. |
| [Grouped-factor research](codex/GROUPED_FACTOR_RESEARCH.md) and [angle-gauge research](codex/ANGLE_GAUGE_RESEARCH.md) | Connect modeled shapes to fit/quantization outcomes and explain why lookup batching/encoding changes the complete resource cost. Some checkpoint figures predate the final native result above. |

The discrete-width effect is visible in the expanded output: **R24/B64/C32** has 272 inner temporary qubits at seven bits versus **544** at eight bits, while preparation-only Toffolis fall from 7,088 to 6,124. These numbers exclude submitted ANGLE circuits and do not establish physical feasibility.

**Remaining gap:** the located artifacts provide partial scans, not a complete marked Pareto frontier over `(R,B,C,brot,J)` using the current oracle architecture. They do not enumerate independently tunable QROAM exponents. The actual host verifier chooses those internally and exposes no submission knob for `optimize_qubits` or k2; overriding the model locally would not change official scoring. A new frontier scan should retain the host's choices, expose their threshold transitions, and state its ANGLE-cost assumptions before using it to prioritize expensive fits. Accuracy and effective normalization would still need separate evaluation.

### What did not work

| Attempt | Recorded result / reason it was excluded |
| --- | --- |
| Cubic inverse approximation | A corrected circuit passed independent quantum checks, but 94 qubits / 2,458 Toffolis cost more than the selected quadratic implementation. Earlier accumulator overflow was corrected. |
| Goldschmidt correction | Numerical models and reversible square components passed, but no complete competitive circuit was established; products, lookup replay, and phase handling remained costly. |
| Boolean phase synthesis / pebbling / Walsh expansion | The 795-AND model matched the integer reference, but schedules required excessive storage or recomputation; ABC increased the AND count. Walsh accuracy required 2,048 terms and too many rotations. |
| 87-Toffoli DQI kernel | Passed full native verification but scored 1.615697e12, worse than the selected 88-Toffoli kernel because the sampled phase-error penalty changed. This does not prove the smaller kernel intrinsically less accurate. |
| Series, terminal blocks, and six-moment look-ahead | Components passed, but tested complete schedules lost on gates or peak width. An opening look-ahead block cost 28,114T / 970Q versus 20,360T / 943Q. |
| Thin cofactor updates | An arbitrary-input round trip failed at shot 199; preliminary resource counts are invalid as submission evidence. |
| Smaller / grouped chemistry dictionaries | Floating-point improvements did not survive quantization. R28/B64/C32 reached 17.0208 mHa in float, but direct 7-bit quantization reached 391.182 mHa, and covariance refitting still missed the cap. Higher precision lost its lookup-cost advantage. |
| Extra atoms / sparse dictionaries / cheaper field multipliers | Added lookup cost outweighed fit gains, residuals remained invalid, or searches found no multiplier below 42 binary ANDs. Search failures do not prove optimality. |

A late **Wronskian decoder** lead remains open: controlled affine transvections cut a squarefree padding adapter from 28,069 to 7,140 Toffolis per direction and passed component quantum tests. A complete total decoder, phase kernel, and adapters were not implemented. It is research evidence, not a submission improvement. Detailed scope and links are in the [Codex experiment report](codex/EXPERIMENT_REPORT.md) and [continuation notes](codex/CONTINUATION.md).

## Claude: cheaper factor shapes, blocked by quantization

Claude concentrated on chemistry shapes and audited the inherited circuits' resource use. The best package remained **9.0096**.

- **Multi-block chemistry:** two blocks with roughly 351 vectors and 175 columns each promised about a 1.47-fold chemistry gain in the host resource model. Adam fits, exact weight projection, variable projection, rank annealing, eigenvector/ISDF initialization, column clustering, and gradual pruning were tried. Pruning eight vectors per round with leverage scoring was most effective in float: errors reached 12.9 mHa at block size 479 and 14.8 mHa at 447. After 8-bit quantization these became about **25 and 28 mHa**, respectively, failing the 20 mHa gate. Direct splitting was much worse. The local quantizer also performed worse on the inherited basin, limiting conclusions about the architecture itself.
- **Angle precision and smaller single-block fits (estimates):** Claude did not build a 7-bit fit or attempt the smaller single-block shapes. The cost scan predicted about 1.07-fold gate savings at 7 bits, with an expected rounding penalty, and about 823–856 scored qubits at 9 bits versus 526 at 8 bits. Single-block B ≤ 511 promised a 1.38–1.49-fold chemistry gain, but its feasibility was inferred from round-3 accuracy evidence rather than tested locally.
- **Normalization and angle loading:** analysis identified only modest normalization headroom and judged the inherited lookup near its floor. Codex subsequently demonstrated a better signed-gauge/batching circuit, so this was a limit of Claude's tested layout, not a general lower bound.
- **DQI:** profiling located the broad Berlekamp–Massey width peak. Removing stored polynomials required expensive reconstruction; an in-place elimination estimate needed about 57,000 field multiplications per direction. Analysis estimated a roughly 26% failure rate for a qd-style variant; no replacement decoder was built. Lowering the error-weight cap consumed too much failure budget. The cheaper upstream derivative-inverse sweep needed extra symbol storage, losing its advantage under the qubit-squared penalty. A small interpolation-skip saving was analyzed but not implemented.
- **Inverse square root:** narrower workspace did not build in the inherited design; linear/cubic fits, separate term phasing, exponent-controlled output rotations, and reduced rotation precision offered no demonstrated improvement. Codex's later affine storage changes supersede the narrower-workspace conclusion for that inherited layout.

### Additional measured results and evidence limits

Claude's expanded [findings](claude/FINDINGS.md) distinguish measured runs, host verification, and estimates. Almost all compute went into chemistry; DQI and inverse-square-root work reproduced and profiled the inherited circuits.

| Chemistry experiment | Float two-body error | What it establishes |
| --- | ---: | --- |
| L-BFGS, R1/B719/C350 | 21.0 mHa | From-scratch single-block control |
| L-BFGS, R2/B351/C175 | 21.1 mHa | Two blocks matched that control at equal effort, but both failed accuracy |
| Variable projection, R2/B351/C190 | 20.3 mHa | More columns approached the cap before quantization |
| Full-rank variable projection, R2/B351/C351 | 17.2 mHa | A float fit passed the physical error threshold; this was not a quantized submission |
| Relaxed inherited R1 solution | 11.3 mHa | A substantially better basin already existed than the from-scratch fits found |

Slow pruning mattered: at B479, dropping 24 vectors per round left 22.5 mHa, dropping eight with a diagonal/row-norm criterion left 15.4 mHa, and eight with leverage scoring left about 13.0 mHa. Quantization then erased that advantage: annealing alone gave 118–128 mHa at B479/B447, and integer descent recovered only 25.0/27.8 mHa. Gradient-sign, curvature-aware, and pair-move quantizers improved successive trial residuals from roughly 163 to 55 to 47 mHa, but these were different trial fits, not a controlled comparison. The inherited-basin calibration still ended at 22.7 mHa versus round 3's approximately 19.6 mHa, identifying quantizer quality as a separate limitation from the float basin.

The R2/B351/C175 host run established **resource counts only**: its unconverged plumbing fit had about 240 mHa error and failed accuracy. The apparent 1.47-fold gain was never a verified usable chemistry improvement. The constrained normalization refitter likewise only reproduced about 48.889 on the inherited solution.

DQI profiling measured 159,664 forward and 133,839 reverse Berlekamp–Massey Toffolis. It also separated the 80-Toffoli Forney point kernel into an 11-Toffoli zero test, 66-Toffoli inverse-parity computation, and three final parity Toffolis. The reported 1,428 live-qubit peak occurs at round 68, while the built circuit reports 1,429 qubits. Inverse profiling measured stage costs of 210 normalization, 308 lookup, 1,043 Horner, and 228 output-shift Toffolis, plus 34 rotations. Alternatives to those circuits remained estimates.

### Corrected conclusions and process lessons

Claude's new [MISTAKES.md](claude/MISTAKES.md) explicitly retracts the claimed overall ceiling near 9.05, the chemistry lookup floor, the 94-qubit inverse limit, the roughly 1,300-qubit decoder floor, and the claim that normalization was exhausted. Codex's verified results are counterexamples. Claude also withdraws the earlier unsupported suggestion that 9.5–10 was realistic, which influenced the user's revised target.

The recurring error was treating a builder constraint, one layout's cost model, or an unsuccessful search as a general limit. Dismissing small gains also prevented testing changes that could accumulate. The revised handoff recommends starting from Codex's selected inverse and signed-gauge lookup scheme, improving quantization from the inherited chemistry basin, and evaluating compact DQI designs with their full repetition penalty. See the [tool index](claude/README.md) for retained scripts.

## AGY: faster chemistry searches, but no lower verified cost

AGY built an analytical offset optimizer reported to reduce fitting time from about 11 minutes to **4.4 seconds**, a quartic tensor-power atom generator, and an automated compile/verify pipeline. These made exploration faster without improving the packaged **9.0096** score.

| Chemistry candidate | Native verification | Cost | Interpretation |
| --- | --- | ---: | --- |
| Incumbent B719/C350, 8-bit | Pass | **757,163,489.50** | Retained |
| C347 | Pass, 50,880 shots | 765,208,795.27 | Fewer gates, worse normalization |
| C345 | Pass, 50,880 shots | 768,464,320.38 | Same trade-off, larger loss |

Reducing columns saved only 0.47–0.78% in the gate-width product while effective normalization rose 1.54–2.30%, worsening the actual score. Both candidates were physically valid; neither was an optimization of the scored objective.

### Additional inverse and decoder investigations

AGY's expanded [experiment report](agy/EXPERIMENT_REPORT.md) describes working-precision sweeps at 28–29 bits, splitting high coefficient bits, early curvature uncomputation, and hierarchical pebble schedules. It reports candidate widths of 86–90 qubits, with a 90-qubit example costing about 2,215 Toffoli-equivalents and 199,350 overall, worse than the inherited 190,712 cost. These are **agent-reported experimental figures**, without a linked candidate-specific native verification result in that section. Its 93-qubit incumbent accounting also differs from the 94-qubit accounting used in the retained round summaries. The common verified cost is unchanged; this report does not treat those width figures as a newly verified milestone or adopt its asserted global polynomial Pareto optimum. Codex's selected 90-qubit circuit demonstrates that smaller width can pay off with a different storage schedule.

For DQI, AGY compares the upstream precomputed derivative-inverse sweep (about 444,030 Toffolis / 1,837 qubits) with the inherited stream decoder (645,755 / 1,429). At the assumed inherited repetition factor, the fast-sweep design needs **fewer than about 1,723 qubits** to break even. This calculation explains why adding storage for a cheaper sweep can lose under the quadratic width penalty; it does not establish a bound on different decoder designs. AGY also reports Fitzpatrick EEA degree-halting trials on 256–1,024-shot sets. No better full decoder emerged: its analysis found that branching/history needed for reversible variable-depth execution could consume the anticipated savings.

### Atom-generation progress and unsuccessful routes

The expanded report records a useful intermediate success: ten tensor-power iterations reduced the quantized reconstruction floor from **19.4272 to 19.1568 mHa** in 8.8 seconds. This was a residual-fitting result, not a fully verified lower-cost submission. It should be distinguished from the unsuccessful replacement and expansion tests below.

Other unsuccessful routes:

- Splitting the dictionary into two disjoint blocks removed needed cross-block correlations; the tested float ALS fit stalled at **479.36 mHa**.
- Dropping five low-weight atoms worsened the reconstruction floor from **19.427 to 19.961 mHa**. Five greedy replacement atoms recovered only 0.077 mHa, reaching **19.885 mHa**, worse than the original floor.
- Expanding from 719 to 725 atoms increased Toffolis from 29,440 to 29,796; small normalization gains could not cover the added preparation cost.

AGY also proposed packing both locator and auxiliary polynomial state into syndrome storage to approach 938 DQI qubits while retaining the fast upstream sweep. Its projected 14–17-fold overall scores are **roadmap estimates**, not verified results. The [scaling analysis](agy/SUBMISSION_AND_SCALING_ANALYSIS.md) includes future-tense verification and optimistic factorization targets; the actual retained score is established by its score summary. See [chemistry conclusions](agy/chemistry_optimization_conclusions.md) for measured results and the [expanded report](agy/EXPERIMENT_REPORT.md) for the broader investigations.

AGY's additional proposed chemistry floor of about 749.2M assumes the inherited gate-width product (29,440 × 526) and an asserted normalization bound. Its claim that going below 740M strictly requires fewer dictionary rows or columns does not apply when the lookup circuit itself changes: Codex reached 708.55M with signed-gauge oracles. Likewise, the failed block-splitting fit is evidence against that tested partition, not proof that every post-hoc partition must fail. The roadmap's overall-score claims also require checking: a 2.5–3.3-fold DQI-only improvement from 9.0096 would yield about **12.23–13.41** overall, not 14.5–17.0, with the other tasks fixed. Larger scores require additional demonstrated gains.

Reusable AGY outputs include the [analytical offset solver](agy/scripts/fast_optimize_offsets.py), [candidate preparation](agy/scripts/prepare_cand.py), [tensor-power generator](agy/scripts/test_power_atoms_on_incumbent.py), and [end-to-end checkpoint verifier](agy/scripts/test_checkpoint.py). Candidate-specific native evidence is retained at [C347 verification](agy/runs/c347_eval/verify.json) and [C345 verification](agy/runs/c345_eval/verify.json).

## Lessons from the round

The successful changes optimized the full scored objective: lookup width and equivalent angle representations mattered as much as arithmetic gate count. Floating-point chemistry fits frequently failed after low-bit quantization. DQI storage reductions paid substantial gate and repetition penalties. Valid candidates could still score worse.

Codex's fixed-shift chemistry analysis bounds effective normalization below by **46.78618**, leaving less than 3.87% improvement from normalization alone relative to its selected model. Additional rank witnesses exclude specific undersized factor families. These are scoped mathematical bounds; they do not prove that other circuits, factorizations, or decoders cannot reach 14.0. Likewise, the other agents' failed fits and architecture estimates should not be treated as impossibility results.

The round ends with a verified incremental improvement and reusable arithmetic/search components, but no demonstrated route to the target. The clearest unfinished lead is the cheaper Wronskian padding component; realizing its benefit still requires a complete reversible decoder and native scoring.
