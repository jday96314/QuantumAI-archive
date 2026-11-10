**Working directory:** `initial_experiments/5/chem-qrom`

## Objective

Reduce the actual scored cost of the **already accepted chemistry factorization** through substantially better synthesis of its ANGLE lookup circuits.

Do not attempt to find a new tensor factorization initially.

This experiment deliberately avoids the difficult 20 mHa numerical-fitting problem by holding the accepted physical representation fixed.

The central question is whether the current angle words admit a significantly cheaper reversible lookup architecture.

## Starting point

The selected chemistry implementation is:

`initial_experiments/4/codex/candidates/norm_cycle2_signed`

Its verified parameters are:

- R=1, B=719, C=350
- Eight-bit angles
- Batch size J=7
- Eleven ANGLE oracles
- Lambda_eff=48.594175
- T_total=30,504
- Q_total=478
- Chemistry cost=708,547,388.96

The Toffoli breakdown is:

- ANGLE lookups, forward and reverse: 10,604
- Inner state preparation: 17,333
- Outer state preparation: 317
- Givens rotations: 2,250

The ANGLE oracles use 123 qubits each, within the 128-qubit cap.

The scored chemistry width includes 141 persistent qubits and a maximum temporary workspace of 337 qubits.

Inner preparation currently sets the temporary-width maximum.

**Consequently, reducing ANGLE scratch alone will not necessarily reduce Q_total.**

Read the accepted round-4 signed-gauge implementation, `chem-qfit` oracle builders, and the completed `inv-dyadic` report before beginning.

The inverse-square-root agent developed a streaming ANF-erasure component that reduced temporary lookup workspace without increasing Toffoli count. Investigate whether that component enables different chemistry lookup strategies.

## Research direction A: Content-aware Boolean synthesis

The current circuits already exploit signed-angle gauges and correlated payload bits.

Investigate whether additional structure exists in the actual table data.

Possible techniques:

- Shared multi-output Boolean subexpressions.
- Reduced decision diagrams.
- Repeated-subtree elimination.
- Content-pruned lookup tries.
- Affine or low-degree Boolean preprocessing.
- Shared prefix computations across output bits.
- Compact XOR networks combined with sparse nonlinear correction terms.
- Measurement-based erasure with deferred phase corrections.

Build the actual resulting circuits and measure compiled Toffolis and qubits.

Do not rely exclusively on abstract Boolean gate counts.

The selected implementation is already highly optimized, so proposed improvements need convincing complete-circuit measurements.

## Research direction B: Address permutation

The dictionary rows may be permuted if the corresponding factorization weights and angle-address mapping are transformed consistently.

Investigate permutations that make nearby table addresses share more bit patterns.

Possible ordering heuristics include:

- Greedy Hamming-distance minimization.
- Trie/subtree similarity.
- Correlated angle-word grouping.
- Clustering by shared Boolean features.
- Ordering optimized directly for the chosen lookup compiler.

One-body address semantics must remain correct.

For two-body rows, update U and W consistently so that the represented Hamiltonian remains unchanged.

Compare the actual compiled cost of the permuted lookup family, rather than assuming that lower Hamming distance produces fewer Toffolis.

## Research direction C: Equivalent angle encodings

The round-4 signed-angle-gauge transformation made two output bits dependent and enabled larger lookup subtrees.

Systematically search other mathematically equivalent angle encodings.

Potential transformations include:

- Sign-equivalent vector representations.
- Equivalent Givens-angle gauges.
- Address-dependent gauge choices where allowed by the representation.
- Affine transformations of stored angle words followed by inexpensive reversible decoding.

Verify that every transformed word sequence decodes to the intended physical projector.

Do not change the physical Hamiltonian merely to simplify the lookup.

The goal is to make the table more compressible while preserving the existing accepted factorization.

## Research direction D: Streaming erasure and batch size

Investigate whether the streaming ANF-erasure technique permits more efficient ANGLE batching.

Explore different J values and corresponding oracle families.

Larger batches may reduce the number of oracles but increase persistent word storage and Q_total.

Smaller batches may reduce persistent storage but increase lookup passes.

Optimize the complete host cost, including:

- ANGLE Toffolis.
- Persistent angle-word storage.
- Maximum temporary workspace.
- Host state-preparation costs.
- Givens costs.
- Oracle count and batch size.

Do not treat the 128-qubit per-oracle cap as the complete chemistry-width constraint.

## Required verification

Every oracle family must implement the exact XOR-difference interface.

Oracle 0 XORs in the first word batch. Subsequent oracles XOR the difference between consecutive batches.

The reverse family must restore the incoming word registers.

Test:

- Every valid address.
- Arbitrary incoming target words.
- Forward and reverse execution.
- Clean scratch.
- No unwanted phase.
- Unchanged selection registers.
- Correct padding and signed-word semantics.

No individual oracle may exceed 128 qubits.

For promising candidates, run the unchanged full chemistry verifier, including its native 64-shot-per-address checks.

If dictionary rows are reordered, verify that the submitted NPZ and regenerated angle words preserve the physical reconstruction and normalization.

## Performance objectives

At unchanged Q and normalization:

- A 1.10× chemistry improvement requires approximately 2,773 fewer total Toffolis.
- A 1.20× improvement requires approximately 5,084 fewer total Toffolis.

Since the current ANGLE family costs 10,604 Toffolis, these correspond to approximately 26% and 48% ANGLE reductions respectively.

These are ambitious targets, not estimates of attainable compression.

Even a smaller verified improvement is useful.

If changed batching also reduces Q_total, account for that directly rather than using the fixed-width thresholds above.

## Research persistence

Explore at least two substantially different lookup synthesis strategies.

Do not declare the table incompressible based solely on one trie or Boolean compiler.

However, avoid unlimited generic Boolean optimization when the measured compiled Toffoli savings are negligible.

Preserve working compiler improvements even if the complete chemistry score does not improve.

## Deliverables

Maintain `RESULTS.md`.

Provide:

- Baseline oracle reproduction.
- Content/statistical analysis of the actual angle words.
- New compilation strategies.
- Exact gate and width comparisons.
- Row/gauge transformations and their inverses.
- Native oracle correctness tests.
- Full chemistry verification results for competitive candidates.

Write only inside `initial_experiments/5/chem-qrom`.

Completed prior work is read-only. Do not inspect concurrently active agent directories.

Do not modify the host verifier, accepted source tensors, or scoring rules.

**The objective is a cheaper verified chemistry circuit without requiring a new physical factorization.**

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