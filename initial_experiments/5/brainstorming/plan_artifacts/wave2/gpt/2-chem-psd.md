# Round 5 Continued — Block-PSD Chemistry Factorization

You own `initial_experiments/5/chem-psd`.

Investigate fundamentally different chemistry factorization architectures, emphasizing **block-positive-semidefinite optimization, increased column capacity, and overlapping nonorthogonal dictionaries**.

The objective is to find a physically accurate factorization with substantially lower host-scored preparation resources.

Unlike the previous multi-block experiments, do not rely primarily on joint approximate diagonalization or short alternating-gradient fits.

## Background

The verified round-4 incumbent has:

- R=1, B=719, C=350
- Eight-bit Givens angles
- T=30,504
- Q=478
- Effective normalization=48.594175
- Chemistry cost=708,547,388.96
- Two-body error=19.9953 mHa

The previous `chem-blocks` and `chem-blocks/jad_search` experiments explored shapes including:

- R4/B180/C88
- R4/B224/C88
- R2/B351/C175

They identified attractive conditional resource costs but generally failed to reach the required 20 mHa two-body accuracy.

For example, the best tested R4/B180/C88 continuous fits remained around 55 mHa.

Those results are credible failures of the tested numerical architectures. They do not establish that other block factorizations are infeasible.

In particular, increasing the number of columns per block and changing the dictionary geometry have not been systematically explored.

Read both prior multi-block reports, their resource-frontier implementations, the completed `chem-qfit` report, and the host physics/scoring code.

## Primary research idea: Optimize block covariances directly

Let P_r contain the symmetric rank-one projector features associated with the vectors in ensemble r.

Each block contributes:

H_r = P_r K_r P_r.T

where:

K_r = W_r W_r.T

and therefore:

K_r >= 0, rank(K_r) <= C.

For fixed dictionaries P_r, first solve the relaxed problem:

minimize ||H_target - sum_r(P_r K_r P_r.T)||_F²

subject to:

K_r >= 0.

**This is a convex optimization problem when the dictionaries are fixed and the rank constraints are removed.**

Implement an efficient projected-gradient, accelerated proximal, operator-splitting, or other appropriate solver.

Use the host's actual symmetric-pair feature metric, including sqrt(2) off-diagonal scaling, the full signed target spectrum, and native reference tensor semantics.

Avoid constructing unnecessarily large dense intermediates where matrix-free contractions or cached Gram matrices suffice.

Verify objective calculations against direct tensor reconstruction.

Track optimization residuals and optimality diagnostics. Do not call a numerical plateau a mathematical lower bound without a sufficiently convincing certificate.

## Experiment A: Determine whether column rank is the bottleneck

For several promising resource-favorable dictionaries, compare:

1. Fully relaxed block-PSD solutions without rank constraints.
2. Rank-88 factorizations.
3. Rank-128 factorizations.
4. Rank-160 factorizations.
5. Higher ranks if the host resources remain competitive.

Candidate families include:

- R4/B180
- R4/B224
- R2/B351

Use the exact host resource model to determine which column counts are worth investigating.

The previous resource scans provide a starting point, not definitive physical feasibility.

If a fixed dictionary's full-rank block-PSD optimum remains substantially above 20 mHa, increasing C cannot solve its problem. That is evidence to change the dictionary, not an invitation to run more weight optimization.

If the relaxed optimum passes but the rank-constrained fit fails, determine the minimum viable C and its complete host cost.

## Experiment B: Overlapping nonorthogonal dictionaries

If rank expansion is insufficient, change the dictionary geometry.

Do not require each group of generator matrices to share an approximate eigenbasis.

Instead investigate dictionaries containing:

- High-leverage vectors shared across ensembles.
- Ensemble-specific vectors derived from residual modes.
- Vectors from the accepted incumbent.
- Nonorthogonal and overcomplete projector bases.
- Newly optimized unit vectors selected to reduce the full physical residual.

Allow controlled redistribution of dictionary capacity between ensembles, subject to the verifier's uniform-shape requirements.

Test at least two materially different initialization strategies.

Use continuous optimization to improve candidate geometry, followed by exact eight-bit quantization and refitting.

Do not infer that a floating-point fit will remain accurate after rounding.

## Experiment C: Physical accuracy and normalization

For promising candidate geometries:

- Restore rank-C block factorizations.
- Optimize one-body identity offsets.
- Optimize the actual SOS gap and effective normalization.
- Quantize the vectors using integer Givens words.
- Refit the weights and offsets against the decoded vectors.

Preserve the verifier's fixed symmetry shifts.

Use the true two-body residual and all physical constraints.

Acceptable numerical residuals alone are insufficient: a candidate must also achieve a competitive lambda_eff×T×Q product.

## Resource and accuracy objectives

The initial target is a complete chemistry improvement of at least 1.3×.

Stretch target: 1.5–1.7×.

The following previously measured resource-only projections illustrate potentially interesting regions:

- R2/B351/C175: approximately 1.45× at incumbent normalization.
- R4/B180/C88: approximately 1.79× at eight bits.
- R4/B224/C88: approximately 1.54× at eight bits.

These are not physical-fit predictions.

Increasing C will raise resource costs, so determine whether the added rank capacity justifies the expense.

## Research expectations

Do not stop merely because an initial JAD or Adam fit exceeds 20 mHa.

The main purpose of this experiment is to distinguish:

1. Dictionary-span limitations.
2. Insufficient block covariance rank.
3. Numerical optimization failures.
4. Quantization failures.
5. Normalization/resource failures.

If the first family is unsuccessful, investigate another resource-favorable family with different dictionary geometry.

An accurate numerical certificate showing why a specific representation fails is useful, but does not establish a global chemistry limit.

## Deliverables

Maintain `RESULTS.md` and preserve:

- Exact host resource scans.
- Fixed-dictionary PSD solver and convergence evidence.
- Rank-versus-error curves.
- Dictionary construction methods.
- Continuous and quantized checkpoints.
- Exact physical reconstruction metrics.
- Full ANGLE circuits and native verification for viable candidates.

Write only inside `initial_experiments/5/chem-psd`.

Read completed round-5 experiments freely, but do not inspect concurrently active agent directories.

Do not modify `from_host/`, previous candidates, or official scoring rules.

**The central objective is to discover whether different block covariances and overlapping dictionaries can make a substantially cheaper chemistry representation physically viable.**