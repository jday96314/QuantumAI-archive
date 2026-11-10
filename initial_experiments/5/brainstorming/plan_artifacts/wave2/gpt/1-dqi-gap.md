# Round 5 Continued — Gap-Aware Compact DQI Decoder

You own `initial_experiments/5/dqi-gap`.

Your objective is to finish the promising DQI exception-repair work from `initial_experiments/5/dqi-repair`, then investigate whether a unified gap-aware decoder can improve both correctness and arithmetic cost.

This is an implementation and optimization task, not merely a feasibility assessment. Continue exploring alternative implementations when an initial design fails, rather than stopping after the first unsuccessful integration attempt.

## Background

The round-4 incumbent has:

- Q = 961
- Mean Toffolis = 1,019,865.8565
- Repetition multiplier R = 1.696846
- 380 phase failures / 102,400 shots
- DQI cost = 1,598,204,229,385

The previous `dqi-repair` agent identified the cause of essentially all these failures: unexpected degree gaps in the compact Euclidean decoder.

Its independent classical model reproduces exactly the original 380 phase failures.

Its repaired classical decoder reduces those failures to **2**, including corrections for immediate completion and terminal degree gaps.

It also implemented and verified a reversible degree-(2,1) correction block. A full forced-stage diagnostic correctly handled 448 targeted exception shots and required 10,506 fewer static Toffolis than the corresponding incumbent path.

However, the forced correction fails ordinary inputs because the agent never implemented general quantum-controlled dispatch.

The previous work is therefore incomplete, not disproven.

Read its `RESULTS.md`, implementation files, classical census, reversible components, and native verification logs before beginning.

## Primary objective: Complete general exception dispatch

Construct a complete reversible decoder that handles:

1. Ordinary degree-one quotient updates.
2. Isolated degree-two quotient updates at arbitrary stages.
3. Immediate completion following an exceptional quotient.
4. Degree gaps near the final two-slot boundary.
5. Repeated degree gaps, to the extent economically justified.

The last category may retain a small failure probability if correcting it would cost more than the corresponding repetition penalty.

The repaired decoder must preserve the original syndrome registers, clear scratch, and apply the correct phase under arbitrary classical measurement masks.

### Investigate two architectures

**Architecture A: Quantum-selected branch routing**

Integrate the existing regular and exceptional blocks using reversible flags, controlled permutations, and shared arithmetic.

Avoid duplicating large arithmetic networks for ordinary and exceptional cases.

Investigate reconstructing branch information from surviving polynomial invariants during uncomputation.

**Architecture B: Unified gap-aware three-slot update**

Instead of dispatching between two different updates, derive a common three-slot transformation that handles ordinary and degree-two quotient cases.

Look-ahead Padé / Euclidean methods may provide useful mathematical inspiration.

The objective is to share field multiplications and cofactor updates between the two cases, using inexpensive conditional routing to select the appropriate data.

This architecture is particularly interesting because the existing exceptional arithmetic block was already cheaper than its ordinary counterpart in one forced-stage experiment.

Do not assume that a rarely activated quantum-controlled branch has negligible gate cost. The host charges its arithmetic regardless of how frequently the quantum controls are satisfied.

## Implementation procedure

First reproduce the existing classical diagnosis and native correction-component tests.

Then build a unified normal/exception component and test it independently before integrating the full decoder.

Test exceptional gaps at every relevant stage, rather than just stage 10.

Include zero-syndrome, one-error, weight-35, immediate-completion, and ordinary random-input controls.

Use the existing million-case classical census to estimate the residual decoding-error distribution before expensive native verification.

Once component correctness is established, integrate the complete forward decoder, phase extraction, and reverse decoder.

Measure the actual live-qubit peak and mean executed Toffolis.

Finally, run the unchanged native DQI verifier for all 102,400 shots, including independent validation salts.

## Resource objectives

If the repaired decoder reaches R approximately 1.206:

- At Q=961, T<=1.148M gives approximately 1.25× DQI improvement.
- At Q=980, T<=1.104M gives approximately 1.25× improvement.
- At Q=1000, T<=1.060M gives approximately 1.25× improvement.

These are conditional budgets, not guaranteed results.

Prefer a design yielding at least 1.25× improvement.

A 1.35–1.40× improvement is a reasonable stretch objective if dispatch can be implemented with little additional workspace.

Do not impose an arbitrary hard width limit. Optimize the actual R×Q²×T objective.

## If the first implementation fails

Do not stop immediately.

Identify whether the obstacle is:

- Insufficient reversible information.
- Incorrect branch cleanup.
- Excessive conditional arithmetic.
- Live-register allocation.
- Exceptional-case handling.
- A genuine mathematical incompatibility.

Then investigate at least one substantially different routing or unified-update design.

If neither is competitive, preserve the working components and provide an explicit resource lower bound or detailed architectural obstruction, clearly distinguishing a proven bound from an implementation limitation.

## Deliverables

Maintain `RESULTS.md` with experiments, failures, resource counts, and reproduction commands.

Preserve:

- Classical reference implementations.
- Reversible component circuits.
- Full candidate builders.
- Register-lifetime profiles.
- Independent verification results.
- The best complete candidate, if one improves the incumbent.

Do not modify the host verifier or previous candidates.

All new files must remain inside `initial_experiments/5/dqi-gap`.

Completed experiments in earlier rounds and round 5 are read-only references. Do not read or modify other concurrently active agent directories.

Use the actual host scoring function and confidence calculation. Never select a candidate based on a favorable circuit-hash sample.

**The objective is to complete a genuinely competitive gap-aware decoder, not merely demonstrate that the previously implemented exception block works.**