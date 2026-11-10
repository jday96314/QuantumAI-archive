# Optimization iteration 4

Three agent workspaces explored the fourth round. **Codex improved the selected local score from 9.0096229084 to 9.3677221531 (+3.97%)**. Claude and AGY retained the inherited package after their alternatives failed accuracy or scored worse. The 14.0 target was not achieved; Claude's subsequent 9.4 target was also not achieved.

The consolidated [round-4 report](../initial_experiments/4/REPORT.md) contains the cross-agent comparison, detailed measurements, corrections, and links. This page summarizes the ideas and their practical outcomes.

| Task | Inherited cost | Best final cost | Cost reduction |
| --- | ---: | ---: | ---: |
| Inverse square root | 190,711.9555 | 182,956.5272 | 4.07% |
| DQI | 1.612737419463 trillion | 1.598204229385 trillion | 0.90% |
| Chemistry | 757,163,489.50 | 708,547,388.96 | 6.42% |

## Codex: bank smaller gains across all tasks

**Inverse square root:** subtracting an affine function of the interval address compressed the lookup constant from 27 to 23 bits. Delayed coefficient loading, early curvature release, shorter offset restoration, and measurement-based output-shift cleanup reduced width from 94 to **90 qubits**. The circuit passed **200,540,928 native shots**, using about 1,782 mean Toffolis and 34 rotations. Cubic interpolation was correct after fixing overflow but more expensive; Goldschmidt, Boolean synthesis/pebbling, and Walsh models did not yield a complete competitive circuit.

**DQI:** blocked compact Euclidean updates consumed syndromes into polynomial cofactors, switched to serial updates late, and deferred completion handling. A cubic GF(4096)-over-GF(16) trace kernel reduced the phase scan from 92 to 88 Toffolis per point. Width fell from 1,429 to **961**, but gates rose to about 1.020 million and repetition to **1.69685**. Full verification had **zero value errors and 380 phase errors in 102,400 shots**. Most squared-width savings were consumed by these increases, leaving only 0.90% cost improvement.

A full 87-Toffoli kernel variant passed but scored worse on its changed sample penalty. Series/terminal schedules and look-ahead components passed component checks but did not improve total resources. A thin-update round trip failed and its preliminary counts were rejected. Field-multiplier searches found no sub-42-AND replacement; this is a search outcome, not an optimality proof.

**Chemistry:** signed angle gauges made two payload bits dependent. Eleven seven-word oracle batches then exploited larger lookup subtrees under the 128-qubit cap. Scored width fell from 526 to **478**. Integer dictionary and weight/offset refits reduced effective normalization to **48.594175**; despite rising to 30,504 Toffolis, total cost fell 6.42%. Full verification passed **50,880 shots**, with physical errors **11.9863/19.9953 mHa**. Supplemental incoming-target checks passed for all eleven oracles, and independent rebuilds reproduced selected outputs.

Smaller/grouped dictionaries again failed quantization: a 28-group float fit reached 17.0208 mHa, but seven-bit rounding gave 391.182 mHa. Extra atoms or precision often lost their resource advantage. A normalization lower bound of 46.78618 leaves limited normalization-only headroom, but does not bound all alternative architectures.

## Claude: better float fits did not survive the local quantizer

Most compute went into two-block chemistry. From-scratch single-/two-block L-BFGS controls reached about **21.0/21.1 mHa**, while relaxing the inherited solution reached **11.3 mHa**. Warm-started pruning was more effective: leverage-based eight-vector steps reached about 13.0 mHa at block size 479, compared with 22.5 mHa for aggressive 24-vector drops. Yet the B479/B447 fits ended at about **25.0/27.8 mHa** after quantization, failing the 20 mHa gate.

The local quantizer also degraded an inherited-basin control to 22.7 mHa versus round 3's approximately 19.6, so these failures implicate the quantizer and basin rather than proving blocked shapes impossible. A two-block pipeline confirmed attractive resource counts only on an inaccurate plumbing fit. Seven-bit and smaller-single-block conclusions were estimates, not newly built fits.

DQI and inverse work reproduced and profiled the existing designs without constructing replacements. The expanded [findings](../initial_experiments/4/claude/FINDINGS.md) and [correction record](../initial_experiments/4/claude/MISTAKES.md) retract unsupported general limits on oracle cost, inverse width, compact decoding, and the overall ceiling near 9.05. Codex's selected circuits are counterexamples. The process lesson is to distinguish a constraint of one builder from a mathematical bound and to test small gains that can accumulate.

## AGY: faster searches and valid but more expensive fits

An analytical SOS-offset gradient solver reportedly reduced fitting from about 11 minutes to **4.4 seconds**. A tensor-power generator improved a reconstruction floor from **19.4272 to 19.1568 mHa** in 8.8 seconds, but no lower-cost packaged candidate resulted. Dropping and greedily replacing five atoms ended at 19.885 mHa, worse than the original floor; adding atoms increased preparation costs.

Column truncations **C347/C345** passed all 50,880 host shots, but scored **765.209M/768.464M**, worse than 757.163M. Gate-width savings of 0.47–0.78% were outweighed by normalization increases of 1.54–2.30%. A tested disjoint-block fit stalled at 479.36 mHa.

AGY also reports narrower inverse/pebble candidates that lost on gate cost, and small EEA degree-halting investigations that did not become better complete decoders. The [expanded report](../initial_experiments/4/agy/EXPERIMENT_REPORT.md) does not link candidate-specific native verification for every such figure, and its incumbent width accounting differs from the retained summaries; those results remain attributed reports. Its fast derivative-inverse DQI sweep analysis gives a break-even width near 1,723 under the assumed gate/repetition counts.

Projected 938–1,312-qubit rational-reconstruction architectures and 14–17 scores are **unbuilt roadmap estimates**. Their overall-score arithmetic also needs care: a 2.5–3.3-fold DQI-only gain yields about **12.23–13.41**, with other tasks fixed. AGY's asserted chemistry floor assumes an inherited gate-width product; Codex improved that product without shrinking the physical dictionary.

## Unfinished lead and evidence

A late Wronskian padding component reduced adapter work from 28,069 to **7,140 Toffolis per direction** and passed component quantum tests. A total reversible decoder, compatible phase kernel, and adapters remained unbuilt. This reopens a research route rather than demonstrating a submission gain.

The [Codex experiment report](../initial_experiments/4/codex/EXPERIMENT_REPORT.md) records full verification, source/rebuild evidence, final archive audits, and rejected components. The [consolidated report](../initial_experiments/4/REPORT.md) preserves all three agents' evidence distinctions. Successful incremental gains do not establish global optimality; unsuccessful fits, resource estimates, and bounded searches do not prove the target impossible.
