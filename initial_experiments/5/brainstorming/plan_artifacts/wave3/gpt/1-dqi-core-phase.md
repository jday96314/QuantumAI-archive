**Working directory:** `initial_experiments/5/dqi-core-phase`

## Objective

Develop a fundamentally different reversible Reed–Solomon decoding architecture that **produces a representation suitable for inexpensive phase extraction directly**, rather than decoding into conventional locator/cofactor polynomials and subsequently performing expensive polynomial conversions.

The ultimate goal is to substantially reduce the complete DQI cost:

`R(epsilon_plus) * Q**2 * mean_Toffolis`

This is a mathematical algorithm-design task followed by circuit implementation. Do not limit the investigation to incremental optimizations of the existing compact Euclidean decoder.

## Starting point

The best verified DQI implementation is:

`initial_experiments/5/dqi-gap/candidates/gap_best_n70_m4095.kmb`

Its resources are:

- Q = 1004
- Mean Toffolis = 1,195,802
- R = 1.19126046
- DQI cost = 1,435,930,717,450
- Zero observed failures in three independent 102,400-shot runs

Its static Toffoli breakdown is approximately:

- Decoder core: 824,740
- Phase extraction and associated preparation: 371,150
- Total: 1,195,890

The previous `dqi-phase` experiment also developed independently verified 11-Toffoli and 23-Toffoli pointwise phase kernels.

However, preparing the polynomial representations needed by these kernels was too expensive.

Its half-degree reciprocal-derivative construction required 1,383 qubits. Its modular F-transformation implementations needed even more workspace.

The `dqi-hybrid` experiment established that merely transplanting the cheaper round-3 Berlekamp–Massey arithmetic into the compact decoder is insufficient: the original BM representation requires additional polynomial banks that defeat the Q² objective.

**Read the completed `dqi-gap`, `dqi-phase`, and `dqi-hybrid` reports and their relevant implementations before beginning.**

These results establish limitations of specific representations, not global lower bounds.

## Primary mathematical idea

The compact decoder maintains cofactor polynomials A and B.

For an error at a root β of B, its magnitude can be expressed as:

    e(β) = β^70 / (A(β) * B'(β))

A sufficiently convenient representation would be:

    F(x) = x^70 * inverse(A(x) * B'(x)) mod B(x)

At each root of B, evaluating F yields the error magnitude.

The existing phase agent verified that, when F is available, the required phase can be applied with an approximately 11-Toffoli point kernel.

The problem is obtaining F reversibly without storing enormous additional polynomial workspaces.

**Investigate whether the decoding recurrence itself can maintain F, or another representation from which the final phase is equally inexpensive.**

Do not assume F must be explicitly materialized. Another representation with a cheap phase functional is acceptable.

## Research direction A: Phase-aware Euclidean recurrence

Investigate whether the quotient updates of the compact Euclidean decoder can simultaneously update a suitable phase/evaluator polynomial.

Potential techniques include:

- Maintaining an inverse cofactor incrementally.
- Updating a modular reciprocal alongside the locator.
- Wronskian or Bézout identities.
- Partial-inverse polynomial representations.
- Reversible approximant-basis updates.
- Transposed or adjoint polynomial algorithms.

Look for transformations that maintain the required information within the existing syndrome/cofactor storage.

The key question is not merely whether a recurrence exists, but whether it can be implemented with inexpensive reversible updates and limited additional storage.

Derive its invariants explicitly.

Implement an independent classical reference and compare against exact Reed–Solomon decoding.

Handle generic degree gaps and squarefree-locator assumptions correctly. An approximate decoder may leave extremely rare exceptional paths uncorrected if the complete confidence-adjusted cost improves.

## Research direction B: Different decoding state representations

Investigate alternatives to both:

1. Conventional Berlekamp–Massey with separate locator/template storage.
2. The current compact Euclidean representation with expensive cofactor updates.

A particularly interesting possibility is a representation retaining the information needed for phase extraction in the syndrome registers themselves.

Investigate whether partial-inverse or rational-reconstruction algorithms can update this representation without materializing both complete polynomial banks.

The existing round-3 decoder used approximately 293,503 Toffolis for its forward/reverse core, substantially fewer than the new gap-aware decoder's 824,740.

However, its workspace was much larger.

The objective is to discover a genuinely different arithmetic/storage trade-off rather than combining incompatible components unchanged.

## Research direction C: Direct phase functional

The verifier only requires the parity phase:

    (-1) ** sum_j(popcount(e_j & mask_j))

It does not require outputting the decoded error vector.

Investigate whether the phase can be obtained from syndromes, locator data, and classical masks without constructing a conventional evaluator polynomial.

Previous work demonstrated that classical-mask interpolation can eliminate approximately 70 evaluation positions with no Toffoli cost.

Consider whether a more substantial algebraic transformation can reduce phase-extraction work.

Any classical preprocessing must depend only on legitimately available classical mask data, not on hidden quantum syndromes or the verifier's true error vector.

The mask registers vary per shot and cannot be treated as fixed compilation-time constants.

## Implementation strategy

Start with small independent Reed–Solomon instances to test new mathematical representations.

For each proposed recurrence:

1. Prove or experimentally establish its algebraic invariant.
2. Compare its outputs with an independent exact decoder.
3. Implement a reversible small-instance circuit.
4. Measure forward and reverse arithmetic.
5. Profile the actual live qubit count.
6. Construct a realistic full-size resource estimate.
7. Integrate the full circuit only when that estimate is promising.

Do not rely on operation-count asymptotics or abstract polynomial invertibility alone.

All arithmetic and restoration costs must be included.

In particular, a mathematically invertible transformation is not automatically an inexpensive reversible circuit.

Test the actual syndrome and runtime-mask contracts, including restoration, scratch cleanup, and phase.

## Performance targets

The current gap-aware decoder costs 1.436 trillion.

At Q=1004 and R=1.19126046:

- Approximately 1,087,093 Toffolis gives a 1.10× DQI improvement.
- Approximately 996,502 Toffolis gives a 1.20× improvement.
- Approximately 797,201 Toffolis gives a 1.50× improvement.
- Approximately 597,901 Toffolis gives a 2.00× improvement.

These are conditional targets.

The long-term architectural goal is approximately 400,000–600,000 Toffolis without a substantial increase in qubits.

A solution with somewhat more qubits or a small nonzero error rate is acceptable if the complete R×Q²×T cost is better.

## Research persistence

Do not stop after demonstrating that one particular polynomial transformation is too expensive.

Explore at least two materially different state representations if the first is unsuccessful.

If no complete improvement emerges, deliver the strongest mathematically verified recurrence, reversible component, and credible whole-circuit resource estimate.

Clearly distinguish:

- A proved algebraic obstruction.
- A reversible-storage obstruction for one representation.
- An implementation that is merely too expensive.
- A hypothesis that has not been fully evaluated.

## Deliverables

Maintain `RESULTS.md`.

Preserve all mathematical derivations, reference implementations, reversible components, resource models, and verifier results.

Only claim a task-cost improvement after complete native verification with the unchanged confidence rule and 102,400 shots.

Use independent validation salts.

All new files must remain under `initial_experiments/5/dqi-core-phase`.

Completed prior experiments may be read but not modified. Do not inspect other concurrently active agent directories.

Do not modify the verifier or exploit circuit-hash/sample selection.

**Your objective is a new decoder architecture whose internal representation makes phase extraction inherently cheap.**

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