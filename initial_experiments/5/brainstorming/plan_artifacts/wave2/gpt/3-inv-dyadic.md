# Round 5 Continued — Adaptive Dyadic Inverse-Square-Root Circuit

You own `initial_experiments/5/inv-dyadic`.

Investigate a new inverse-square-root circuit architecture based on adaptive dyadic segmentation, local integer polynomial approximations, and low-workspace reversible arithmetic.

The objective is a substantial reduction in the complete Q×(Toffolis + synthesized_T/4) task cost.

This is an independent architecture experiment, not another general attempt to polish the existing circuit.

## Background

The verified round-4 Codex circuit has approximately:

- Q=90
- Mean Toffolis=1,782
- 34 non-Clifford rotations
- Approximately 2,033 Toffoli-equivalent total work
- Task cost=182,956.5272

The architecture uses normalized input coordinates, piecewise quadratic approximation, carefully optimized fixed-point Horner arithmetic, and measurement-based cleanup.

Previous rounds achieved major breakthroughs by eliminating the original Newton refinement and optimizing coefficient representation and reversible register lifetimes.

However, they have not comprehensively explored adaptive dyadic segmentation as a replacement for the existing interval architecture.

Read the selected round-4 builder, its profiling, the inverse-square-root verifier, and previous approximation experiments.

## Primary hypothesis

The existing interval partition may use more arithmetic precision or lookup work than necessary in some regions.

A dyadic partition could allow:

- Intervals selected directly from binary input prefixes.
- Different interval widths in different parts of the normalized domain.
- Narrower local coordinates.
- Reduced polynomial coefficient ranges.
- Simpler normalization and interval selection.
- More efficient sharing of registers between approximation and phase application.

The critical question is whether these benefits outweigh additional lookup entries and more complicated interval logic.

Do not assume they will.

## Experiment A: Adaptive interval construction

Build an independent exact-integer numerical model.

Partition the normalized mantissa domain into dyadic intervals, with boundaries aligned to binary prefixes.

Generate candidate partitions using curvature/error estimates, then refine them against the verifier's exact relative-error tolerance.

Test at least:

1. Piecewise-linear interpolation with adaptive intervals.
2. Piecewise-quadratic interpolation with adaptive intervals.
3. Mixed-degree approximation, using linear segments where possible and quadratic segments elsewhere.

Optimize polynomial coefficients directly against the integer arithmetic used by the eventual reversible circuit.

Do not rely exclusively on real-valued polynomial approximation error.

Account for coefficient truncation, intermediate rounding, overflow, final phase conversion, and worst-case absolute error margin.

## Experiment B: Reversible arithmetic architecture

For promising approximations, investigate:

- In-place Horner accumulation.
- Borrowing coefficient registers for intermediate products.
- Narrow local-coordinate representations.
- Delayed coefficient loading.
- Early measurement-based coefficient deletion.
- Low-ancilla multiply-add schedules.
- Shared normalization and interval-selection logic.

Measure actual live-register peaks rather than summing independently estimated component widths.

The existing 90-qubit circuit is highly optimized. Avoid introducing a larger accumulator or lookup bank without demonstrating that its arithmetic savings compensate for the width increase.

## Experiment C: Alternative segmentation and arithmetic

If the first dyadic construction is uncompetitive, do not stop immediately.

Investigate at least one materially different variant, such as:

- A prefix tree with asymmetric interval refinement.
- A hybrid of the incumbent centered quadratic and dyadic intervals.
- Reusing input bits as the local coordinate rather than materializing a separate offset.
- A different Horner/carry schedule that prioritizes reducing peak qubits rather than Toffoli count.

Preserve working component improvements even if a complete new segmentation does not win.

## Accuracy and scoring

The verifier checks every x from 1 through 3,133,452, with 64 repetitions to exercise measurement branches.

The relative phase tolerance is 2^-24 / sqrt(x).

A candidate must restore the input, clear scratch, preserve measurement-phase corrections, and satisfy the complete phase requirement.

Use shortened tests and exact integer models for rapid iteration.

Final candidates must pass all 200,540,928 native shots.

Account for rotation synthesis using the verifier's actual worst-case margin; nominal arithmetic accuracy alone does not establish a valid synthesized circuit cost.

## Performance objectives

Current task cost: 182,956.5272.

Initial meaningful target: below 165,000.

Preferred target: below 155,000.

Stretch target: below 140,000.

For example, a hypothetical 80-qubit circuit with 1,900 Toffoli-equivalent work would cost 152,000, an improvement of approximately 1.20×.

This is an illustrative resource combination, not a prediction.

Prioritize measured Q×T-equivalent cost over either Q or gates independently.

## Research expectations

Do not conclude that dyadic interpolation is impossible simply because one particular partition produces too many entries.

Investigate the trade-off between interval count, polynomial degree, coefficient precision, and arithmetic width.

Likewise, if the best numerical approximation produces a costly reversible implementation, determine whether a different register schedule fixes the problem.

Test at least two substantially different architectures before concluding that the explored design family is uncompetitive.

A credible negative result should include an approximation-versus-resource frontier, not merely one failed circuit.

## Deliverables

Maintain `RESULTS.md` with:

- Integer approximation models.
- Dyadic interval definitions.
- Coefficient optimization procedures.
- Worst-case error and synthesis-margin checks.
- Register-lifetime profiles.
- Component and complete circuit costs.
- Native verifier results.
- Reproduction commands.

Preserve the best verified incumbent until a superior complete circuit passes.

Write only inside `initial_experiments/5/inv-dyadic`.

Previous rounds and completed round-5 experiments are read-only. Do not inspect concurrently running agent directories or modify the host verifier.

**The objective is to determine whether adaptive approximation and register scheduling can jointly produce a significantly cheaper inverse-square-root circuit.**