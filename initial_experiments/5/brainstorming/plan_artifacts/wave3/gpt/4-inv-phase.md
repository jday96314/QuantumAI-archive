**Working directory:** `initial_experiments/5/inv-phase`

## Objective

Investigate fundamentally different ways to apply the required inverse-square-root phase **without necessarily materializing a complete fixed-point approximation of 1/sqrt(x)**.

The objective is a substantially cheaper reversible phase circuit.

Do not spend the main effort repeating the adaptive segmentation, integer-coefficient fitting, or Horner register-scheduling searches already performed in rounds 3–5.

## Background

The selected inverse-square-root circuit costs:

- Q=90
- Mean Toffolis approximately 1,782
- 34 non-Clifford rotations
- Approximately 2,033 total Toffoli-equivalent work
- Task cost=182,956.5272

It passes all 200,540,928 native verification shots.

The recent `inv-dyadic` experiment independently generated a 269-leaf binary-prefix approximation that nearly reproduced the incumbent coefficient table.

Its complete implementation also required 90 qubits and essentially the same Toffoli count.

Alternative variable-depth and mixed-linear approximations were more expensive.

The experiment also developed a streaming ANF-erasure component that reduced lookup workspace, but it did not improve the incumbent's total cost.

These results suggest that the existing fixed-point Horner architecture is highly optimized.

They do not establish that the phase itself must be synthesized using that architecture.

Read the selected round-4 inverse builder, the `inv-dyadic` report, previous Goldschmidt/Walsh experiments, and the verifier.

## Primary hypothesis

The task requires the diagonal phase:

    exp(2*pi*i/sqrt(x))

It does not require outputting an approximation of 1/sqrt(x).

Consequently, a circuit may be able to apply the phase incrementally, through a different arithmetic representation, or directly from polynomial and lookup terms.

The goal is to avoid some of the workspace and arithmetic required to construct a conventional complete fixed-point output.

## Research direction A: Incremental phase accumulation

Investigate whether the existing polynomial approximation can be reorganized so phase contributions are applied before the complete result is assembled.

Potential techniques include:

- Applying rotations directly from partial Horner accumulators.
- Computing and phasing individual output-bit groups.
- Digit-serial phase accumulation.
- Reusing arithmetic registers after their contribution to the final phase is applied.
- Measurement-based deletion with correct deferred phase repair.
- Avoiding an independently materialized final phase word.

Any proposed transformation must preserve the exact desired overall phase, including carries, truncation, and intermediate rounding.

A simple decomposition into individually phased polynomial terms is not automatically correct when fixed-point carries or nonlinear products are involved.

Derive the required arithmetic identities explicitly.

Build independent exact-integer models before emitting quantum circuits.

## Research direction B: Phase-oriented approximation

Investigate approximation structures designed directly for phase application.

Examples worth evaluating include:

- A coarse prefix-dependent phase plus a small residual correction.
- Piecewise phase-polynomial synthesis.
- Low-rank or structured decomposition of the phase function.
- Digitwise phase tables combined with inexpensive residual arithmetic.
- Hybrid lookup/arithmetic methods that never construct the entire inverse-square-root value.

Evaluate the complete cost of selecting the approximation, loading any parameters, applying the phase, and cleaning all workspace.

Do not assume that lookup tables or diagonal phase gates are free.

Previous Walsh-based phase synthesis required too many terms. A new approach must explain what structure avoids that failure.

## Research direction C: Alternative arithmetic representation

If direct phase synthesis is initially unsuccessful, investigate a genuinely different reversible arithmetic formulation.

Potential directions include:

- Phase-oriented truncated multiplications.
- Transposed or middle-product arithmetic.
- Bit-sliced accumulation that discards unneeded output information.
- Alternative signed-digit representations.
- Register-borrowing strategies that eliminate a live product bank.
- A phase-aware multiply-add primitive that avoids separately constructing and storing the completed product.

Do not merely rescan the same coefficient widths and interval boundaries.

Focus on reducing the complete arithmetic-and-phase circuit's live state or nonlinear gate count.

## Mathematical and implementation requirements

The verifier's input domain is:

    1 <= x <= 3,133,452

The required relative phase tolerance is:

    2^-24 / sqrt(x)

The complete circuit must restore the input and clear all scratch.

Non-Clifford rotation synthesis must use the verifier's actual accumulated error margin.

The incumbent has relatively little worst-case numerical margin, so a phase decomposition must account carefully for approximation and synthesis errors.

Use an independent exact fixed-point model to test candidates.

For every promising mathematical representation:

1. Derive the complete phase expression.
2. Establish its finite-precision error behavior.
3. Prototype the reversible arithmetic or phase component.
4. Measure Toffolis and peak qubits.
5. Include all forward/reverse cleanup and measurement corrections.
6. Estimate the full task cost.
7. Build and exhaustively verify promising complete candidates.

The final native test must include all 200,540,928 shots.

## Performance objectives

The incumbent cost is 182,956.5272.

Initial target: below 170,000.

Preferred target: below 155,000.

Stretch target: below 140,000.

For example, 78 qubits with 1,800 Toffoli-equivalent work would cost 140,400.

That combination is illustrative, not an expected outcome.

A small component saving is useful only if it survives full width, gate, and synthesis accounting.

## Research persistence

Investigate at least two substantially different phase/arithmetic architectures if the first is uncompetitive.

Do not stop simply because one direct-phase lookup or phase-polynomial decomposition uses too many gates.

However, if the available phase-function structure cannot support meaningful savings, provide explicit resource evidence rather than continuing generic coefficient tuning.

Distinguish failure of a representation from a mathematical lower bound.

## Deliverables

Maintain `RESULTS.md`.

Preserve:

- Mathematical derivations.
- Exact fixed-point phase reference models.
- Approximation/error frontiers.
- Reversible component circuits.
- Live-register profiles.
- Full gate and rotation-synthesis costs.
- Native full-domain verification results.
- Reproduction commands.

Keep the accepted inverse-square-root circuit as the fallback.

Write only inside `initial_experiments/5/inv-phase`.

Completed prior research may be used as read-only reference material. Do not inspect other active agent directories.

Do not modify the host verifier or scoring rules.

**The objective is to discover a different way to implement inverse-square-root phasing, not merely produce another optimized Horner polynomial.**

## Shared Constraints and Research Standards

### 1. Workspace isolation

- Your assigned working directory is the only location where you may create or modify project files.
- Do not read, search, modify, or copy files from other **concurrently active agent directories**.
- Completed experiments from rounds 1–5 may be used as read-only references.
- Treat `from_host/`, the reference tensors, accepted submissions, and competition scoring code as immutable.
- Do not modify shared Git state, install dependencies into shared environments, or create files in another agent's workspace.
- Use experiment-local caches, temporary files, checkpoints, and logs wherever practical.

### 2. Authoritative baselines

The current best verified components are:

| Task | Cost | Source |
|---|---:|---|
| Inverse square root | 182,956.5272 | Round-4 Codex |
| DQI | 1,435,930,717,449.94 | Round-5 `dqi-gap` |
| Chemistry | 708,547,388.9553 | Round-4 Codex |

**Do not accidentally optimize against the superseded round-4 DQI decoder as though it were the current incumbent.**

Preserve the best accepted submission for your task until a complete replacement passes verification and improves its actual scored cost.

### 3. Optimize the real competition objective

Use the unchanged host formulas:

- **Inverse:** `Q * (mean_Toffolis + synthesized_T / 4)`
- **DQI:** `R(epsilon_plus) * Q**2 * mean_Toffolis`
- **Chemistry:** `lambda_eff * T_total * Q_total`

For DQI, use the host's exact one-sided Clopper–Pearson confidence calculation with `delta = 2^-128`.

For chemistry, include all analytically priced state-preparation and Givens costs, actual ANGLE costs, and the host's persistent-plus-maximum-temporary width calculation.

A reduction in gates, qubits, approximation error, or runtime does not necessarily constitute a scored improvement.

### 4. Verification and experimental integrity

- Use abbreviated tests for development, but unchanged full native verification for final claims.
- DQI final verification requires 102,400 shots.
- Inverse-square-root final verification requires the complete 200,540,928-shot workload.
- Chemistry requires complete native physical and ANGLE verification, including its per-oracle 128-qubit limit.
- Use independent validation salts after freezing promising candidates.
- Never select equivalent circuits based on favorable random samples, submission hashes, or verifier salts.
- Never weaken tolerances, modify reference inputs, skip required correctness checks, or change confidence calculations.
- Preserve exact candidate bytes and verification JSONs.

Distinguish explicitly between mathematical correctness, numerical-model success, component verification, estimated complete resources, and a fully native-verified task-cost improvement.

### 5. Research persistence and failure interpretation

**Do not equate an unsuccessful implementation with a fundamental limitation of the underlying algorithm.**

When an approach fails:

1. Identify the specific mathematical, numerical, implementation, or resource bottleneck.
2. Determine whether the failure applies only to the tested representation.
3. Investigate at least one materially different implementation or mathematical formulation when a credible alternative exists.
4. Quantify the resources needed for the approach to become competitive.
5. Preserve working components and mathematical insights even if the complete architecture fails.

A failed optimizer, an insufficient numerical budget, or a local minimum is not proof that a feasible solution does not exist.

Conversely, do not continue refining an architecture after a credible, correctly scoped lower bound rules out competitive resources.

Avoid unsupported claims of global optimality, impossibility, or fundamental precision floors.

### 6. Resource efficiency

Four agents may be operating simultaneously.

- Avoid unnecessarily monopolizing shared CPU, GPU, RAM, or disk resources.
- Prefer inexpensive mathematical and component feasibility checks before constructing enormous circuits.
- Use efficient batched numerical evaluations and existing verified infrastructure where appropriate.
- Monitor resource consumption during large searches.
- Record hardware contention or incomplete numerical convergence when they materially affect experimental conclusions.

Do not sacrifice correctness or research quality merely to finish quickly.

### 7. Progress documentation

Maintain a continuously updated `RESULTS.md` in your assigned directory.

Record:

- Hypotheses investigated and their rationale.
- Important mathematical derivations or invariants.
- Implementations attempted.
- Successful and unsuccessful experimental results.
- Numerical accuracy and actual resource measurements.
- Exact reproduction commands.
- Locations of candidate circuits, checkpoints, and verification JSONs.
- Remaining limitations and the strongest next research direction.

Preserve the best working candidate throughout the experiment.

Before concluding, provide a concise final assessment that separates:

**A. Verified improvement:** A complete candidate passes native verification and beats the incumbent.

**B. Validated research breakthrough:** A mathematical algorithm or reversible component demonstrably improves an important subproblem, but integration remains incomplete.

**C. Scoped negative result:** The investigated implementation is uncompetitive, with concrete evidence explaining why.

**D. Unresolved hypothesis:** The approach remains plausible but was not evaluated sufficiently to establish its merit.

Do not present B, C, or D as A.

### 8. Research objective

The broader objective is to discover algorithmic improvements capable of moving the combined competition score substantially beyond its current level.

The reported public leaderboard leader has reached approximately **18.9**. Treat this as motivation to investigate fundamentally different algorithms, not as evidence for any particular implementation.

Do not assume the current architectures are near globally optimal merely because previous agents exhausted familiar optimization methods.

**Prefer meaningful algorithmic discoveries over cosmetic progress, but bank verified incremental improvements whenever they are available.**