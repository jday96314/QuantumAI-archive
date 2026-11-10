**Working directory:** `initial_experiments/5/chem-tree`

## Objective

Investigate whether substantially different Givens-tree topologies, orbital orderings, and low-precision representations can overcome the quantization barrier that has defeated smaller chemistry factorizations.

This is not another generic floating-point tensor-fitting task.

The primary objective is to change the discrete representation itself, potentially making accurate low-resource factorizations possible.

## Background

The accepted chemistry factorization uses:

- R1/B719/C350
- Eight-bit Givens angles
- A shared 75-edge Givens tree
- Two-body error 19.9953045 mHa
- Chemistry cost 708,547,388.96

Recent experiments found that excellent continuous solutions often become unusable after integer quantization.

Examples:

- R2/B351/C189 has a demonstrated continuous two-body error of approximately 19.979 mHa.
- The rounded eight-bit dictionary derived from that basin has a joint-span lower bound exceeding 140 mHa.
- A saved B479 quantized dictionary has a certified full-PSD residual above 20.156 mHa.
- Other exact-word searches improved numerical errors but often increased normalization dramatically.

These results apply to particular integer dictionaries.

They do not establish that the same factorization shapes are infeasible under different Givens trees or orbital-coordinate representations.

Most prior integer optimization kept the Givens tree fixed.

Read the completed `chem-qfit`, `chem-discrete`, and `chem-psd` reports before beginning.

## Primary hypothesis

A fixed-precision Givens representation defines a discrete set of representable orbital vectors.

Changing the Givens tree changes that discrete set.

A tree that is favorable for one dictionary may be unfavorable for another.

Therefore, the nearest quantized representation under the incumbent tree need not be close to the best representation achievable with a different tree.

The search should optimize the shared tree and the integer dictionary together, rather than treating the tree as fixed infrastructure.

## Research direction A: Tree topology and orbital ordering

First inspect the actual host verifier's tree validation and decoding semantics.

Determine precisely which tree topologies and edge orders are permitted.

Do not assume arbitrary permutations are valid.

Construct several genuinely distinct valid tree families, potentially including:

- Alternative balanced topologies.
- Highly asymmetric trees.
- Orbital-correlation-based trees.
- Amplitude-statistics-based trees.
- Hierarchical orbital clusters.
- Greedy trees optimized for quantization sensitivity.

The tree is shared across all angle-table addresses.

Do not assign a separate tree to every vector.

For each candidate tree, implement correct continuous-angle encoding, integer-word encoding, and physical-vector reconstruction.

Test exact round-trip behavior and independently verify decoded vectors.

Any orbital reordering must be interpreted as a change in representation, not an unauthorized change to the physical reference Hamiltonian.

## Research direction B: Quantization-aware tree scoring

Develop a fast surrogate that predicts which trees are favorable for a collection of physical vectors.

Avoid optimizing only average vector L2 error.

The physical objective depends on rank-one projectors and their contributions to the complete two-body tensor.

Consider sensitivity-weighted quantization errors using:

- Dictionary leverage.
- Factor-column magnitudes.
- The symmetric-pair tensor metric.
- The physical residual's leading modes.
- Jacobians of the decoded projector vectors.
- Nonlinear error cancellation after weight refitting.

For a given continuous dictionary, compare multiple trees using the exact decoded integer vectors followed by rank-constrained weight refits.

Measure actual two-body reconstruction error.

A numerically better tree is useful only if its improvement survives the factorization model.

## Research direction C: Joint tree and discrete-dictionary optimization

Once promising tree families are identified, optimize the integer words directly.

Investigate:

- Tree-edge swaps followed by local angle repair.
- Coordinated subtree restructuring.
- Orbital-order mutations.
- Multi-vector integer angle changes.
- Residual-guided atom replacement.
- Weight and offset refitting after topology changes.
- Annealed or population-based searches across distinct tree families.

Allow the search to move between substantially different discrete basins.

Do not rely only on nearest rounding of a continuous vector.

Preserve the best exact-word dictionary for each tree.

Compare candidates at matched computational budgets.

## Candidate shapes

Use several distinct numerical controls:

1. The accepted R1/B719/C350 dictionary, to calibrate encoding and quantify the effect of changing the tree.
2. At least one smaller single-block shape, such as R1/B639/C350.
3. A two-block continuous basin such as R2/B351/C189.
4. A saved near-valid exact-word two-block candidate, such as R2/B479/C175.

Do not evaluate every tree only against one failing initialization.

For the two-block continuous candidates, determine whether a new tree substantially reduces the extreme quantization error.

A modest improvement in vector representation error may be insufficient; evaluate the physical tensor directly.

## Research direction D: Alternative low-precision representations

If topology changes alone are insufficient, investigate whether the same number of angle bits can be used more effectively.

Possibilities include:

- Different angle gauge conventions.
- Nonuniform angle codebooks, if they can be decoded legally and cheaply.
- Globally coordinated signed encodings.
- Joint optimization of tree topology and quantization grid.
- Representations designed around the physical residual rather than vector reconstruction accuracy.

The host still expects the prescribed integer angle semantics.

Any alternative codebook must therefore be implemented through an allowed representation or an explicitly charged reversible decoding circuit.

Do not silently change the verifier's interpretation of angle words.

## Physical and resource objectives

The native physical requirements include:

- One-body error <=20 mHa.
- Two-body error <=20 mHa.
- Strict Lambda and E_gap constraints.
- Correct fixed symmetry shifts.
- A valid shared Givens tree.
- Correct integer angle words.
- Correct ANGLE oracle behavior and width cap.

Optimize the actual task cost:

    lambda_eff * T_total * Q_total

Do not optimize reconstruction residual alone.

The preferred research result is a new tree/discrete dictionary that makes a resource-favorable smaller factorization physically valid.

A complete 1.20× chemistry improvement would be significant.

A 1.30× or greater improvement is a stretch objective.

A substantial quantization-error reduction on a previously unusable low-cost shape is also valuable research evidence, but is not a verified score improvement.

## Research persistence

If one tree family fails, investigate a genuinely different topology or representation before stopping.

Do not interpret a failed tree optimizer as proof that topology is irrelevant.

Conversely, if several distinct tree families produce essentially identical decoded quantization behavior, investigate whether the failure is caused by deeper limitations of the factorization geometry.

Report that distinction explicitly.

## Deliverables

Maintain `RESULTS.md`.

Preserve:

- Valid tree generators and encoders.
- Independent tree-decoding tests.
- Tree/topology search trajectories.
- Exact integer-word checkpoints.
- Continuous-versus-quantized physical comparisons.
- Weight/offset normalization results.
- Host resource calculations.
- Native verification results for competitive candidates.

Write only inside `initial_experiments/5/chem-tree`.

Use completed prior experiments as read-only references. Do not inspect concurrently active agent directories.

Do not modify the host verifier or reference physical tensors.

**The objective is to escape the current quantization landscape by changing the representation in which the chemistry factorization is optimized.**

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