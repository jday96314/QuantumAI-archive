# Round 5 Continued — Direct Discrete Chemistry Optimization

You own initial_experiments/5/chem-discrete.

Your objective is to discover substantially cheaper, physically valid chemistry
factorizations by optimizing the actual integer Givens-angle representation,
rather than relying on floating-point fitting followed by quantization.

Focus on discrete dictionary geometry, multi-atom changes, normalization-aware
optimization, and escaping local minima.

Do not repeat the previous approach of simply finding better continuous
solutions and rounding them.

## Background

The selected chemistry incumbent has:

R=1, B=719, C=350
brot=8
T=30504
Q=478
lambda_eff=48.594175
cost=708547388.9553
two-body error=19.9953045 mHa

The completed round-5 chem-qfit experiment established several important facts.

1. The accepted integer representation can be reproduced exactly.
2. Floating-point solutions with excellent residuals often become extremely
   inaccurate when converted to integer angle words.
3. Direct integer optimization is much more effective than rounding some
   floating-point basins.
4. Smaller physically valid dictionaries exist, but their normalization
   frequently increases enough to make their scores worse.
5. The strongest smaller two-block integer solutions are close to, but still
   above, the 20 mHa threshold.

Particularly interesting observations:

- R2/B479/C175 floating-point error: approximately 13.15 mHa.
- R2/B479/C175 existing integer solution: approximately 22.62 mHa
  in one direct integer-continuation experiment.
- Other B479 integer basins: approximately 24.95 mHa.
- R2/B447/C175 old quantized dictionary has a relaxed joint-span lower
  bound of approximately 20.596 mHa.
- A B671 row-exchange experiment improved 24.505 to 23.799 mHa.
- The accepted B719 integer basin can be improved to approximately
  19.516 mHa without changing its scored dimensions.

The B447 lower bound applies to that particular dictionary, not to all
integer dictionaries with that shape.

Read the previous chem-qfit RESULTS.md, numerical checkpoints, integer
optimizers, and the complete native scoring implementation.

## Primary hypothesis

A high-quality discrete dictionary may lie far from the nearest quantized
version of the best continuous dictionary.

Consequently, the correct problem is not:

    Find optimal continuous U, then round its Givens angles.

Instead, directly search the space of integer angle words, with the other
factorization variables continuously refitted around each discrete change.

The real optimization objective is the verifier's scored resource product:

    lambda_eff * T_total * Q_total

subject to the exact physical constraints.

## Experiment A: Direct multi-vector integer optimization

Start from the strongest available exact-word solutions, particularly the
near-valid two-block B479 candidates.

Implement a discrete optimizer capable of coordinated changes to several
vectors at once.

Investigate at least two genuinely different methods:

1. Residual-guided multi-atom coordinate search.

   Use physical tensor residuals and their leading modes to propose correlated
   angle-word changes across multiple vectors.

   Refit block weights and identity offsets after promising moves.

   Do not restrict the search to independently improving one angle at a time.

2. Basin-changing integer moves.

   Investigate row exchanges, blockwise replacements, correlated donor/receiver
   changes, integer-angle simulated annealing, or another suitable combinatorial
   search method.

   Permit temporarily worse local moves when there is a credible mechanism
   for escaping the current basin.

   Use exact physical checks to evaluate the resulting candidates.

Preserve the best candidate throughout the search.

Use efficient GPU/batched calculations where possible.

## Experiment B: Quantized atom replacement

The existing chemistry factorization uses rank-one projectors derived from
unit vectors represented by discrete Givens angles.

Generate new candidate vectors directly in that representable space.

Good proposal mechanisms might include:

- Leading residual eigenmodes.
- Quantized combinations of existing dictionary vectors.
- Replacement of correlated or low-leverage atoms.
- Coordinated replacement of multiple atoms.
- Local Givens-tree transformations.
- Residual-driven generation of entirely new integer vectors.

Evaluate candidates using the full symmetric-pair tensor metric.

Refit weights and offsets after replacements.

Avoid relying solely on the improvement in an unconstrained vector-span
proxy: the actual per-block PSD/rank-C structure must also be preserved.

## Experiment C: Optimize normalization simultaneously

The previous B711 and B700 candidates were physically valid but scored worse
because lambda_eff increased significantly.

Do not optimize reconstruction error to 20 mHa and only then discover that
normalization makes the candidate uncompetitive.

Instead, incorporate normalization into candidate selection.

Use a constrained objective that rewards improvements to lambda_eff while
maintaining valid one-body and two-body reconstruction errors.

The solver must respect the exact native fixed symmetry shifts and SOS gap.

Include an explicit safety margin around the strict gap constraint.

Retain physical feasibility throughout the final optimization stage.

## Candidate shapes

Prioritize resource-favorable shapes with credible numerical starting points.

Candidates include:

- R2/B479/C175
- R2/B447/C175
- R1/B639/C350
- Nearby shapes suggested by the existing host resource frontier.

A more aggressive R1/B511 family may be explored if the new discrete optimizer
demonstrates substantial improvements on easier shapes.

For perspective, at the previously measured T/Q resource estimates:

- R2/B479/C175 can beat the current score if lambda_eff is below
  approximately 54.93.
- R2/B447/C175 can beat the current score if lambda_eff is below
  approximately 58.73.
- R1/B639/C350 can beat the current score if lambda_eff is below
  approximately 53.33.

These are conditional break-even thresholds.

Recompute them using actual candidate T/Q values whenever shape, precision,
or lookup construction changes.

A passing physical fit is not sufficient; its actual scored cost must win.

## Experimental process

First reproduce the best saved integer checkpoints from chem-qfit.

Establish matched numerical controls on identical starting candidates.

Run the new discrete optimizers with enough iterations to characterize their
progress, rather than stopping after a short initial plateau.

Compare proposed algorithms at matched computation budgets.

Preserve checkpoint trajectories so a superior intermediate candidate is
not lost if a later optimization becomes unstable.

For every significant checkpoint, record:

- R/B/C and Givens precision.
- Exact quantized two-body error.
- Exact one-body error.
- Lambda and E_gap.
- Effective normalization.
- Estimated full T/Q.
- Estimated task cost.
- Native status, when tested.

Use the actual decoded integer words and host physical metric.

## Verification requirements

Only promote candidates passing:

- Both physical reconstruction thresholds <=20 mHa.
- Strict Lambda and E_gap constraints.
- Exact Givens tree and word encoding.
- Correct XOR-difference ANGLE lookup behavior.
- The 128-qubit per-oracle cap.
- Full forward/reverse native verification.
- A task cost lower than the selected incumbent.

Preserve all intermediate artifacts and complete verification JSONs.

## Performance objectives

Initial target: at least 1.10x chemistry improvement.

Preferred target: at least 1.20x.

Stretch target: 1.35x or more.

Do not stop merely because an initial integer optimizer gets stuck above
20 mHa.

Investigate multiple substantially different move families and initialization
strategies before concluding that the tested approach is ineffective.

If a shape has a convincing fixed-dictionary capacity obstruction, move to
dictionary replacement rather than spending additional time adjusting weights.

Distinguish numerical local minima from mathematical impossibility.

## Deliverables

Maintain RESULTS.md with exact reproduction commands, resource measurements,
successful and unsuccessful hypotheses, and candidate checkpoint locations.

Preserve all improved integer dictionaries, even if they do not yet meet
the physical threshold.

Produce native-verified submission files for any candidate that beats the
current incumbent.

Write only inside initial_experiments/5/chem-discrete.

Completed previous experiments may be read as references.
Do not read or modify directories belonging to other concurrently running
agents.

Do not modify the host verifier, source reference tensors, scoring rules,
or previously accepted submissions.

The objective is to discover a substantially better discrete physical
representation, not merely improve the floating-point fitting pipeline.