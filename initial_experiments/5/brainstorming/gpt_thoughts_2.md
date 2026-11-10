# Reassessment: Claude is too pessimistic about the ceiling, but its resource analysis is useful

**Yes, the additional documents change some of my recommendations. However, I don't think Claude has established anything approaching a 1.6× upper bound for chemistry, much less a ceiling around 16.7 for the overall competition.**

My conclusions are:

1. **Claude is correct about the current chemistry bottleneck.** The detailed resource breakdown makes my earlier suggestion of a relatively straightforward 2× chemistry improvement look too optimistic.
2. **Claude's 1.6× chemistry ceiling is not established.** It's an estimate for a limited family of factorization shapes and lookup implementations. It is not a lower bound on the resources required by the verifier.
3. **The DQI analysis is more actionable than my previous assessment.** Exception repair deserves higher priority, and Claude identifies a particularly interesting way to obtain a division-free phase scan without the enormous additional workspace of the published construction.
4. **I would revise my earlier multi-block chemistry proposal.** Joint diagonalization is still worth investigating, but not as the highest-priority experiment. The host's discrete workspace thresholds make some otherwise appealing shapes substantially less attractive.
5. **I would not accept Claude's 12–14 prediction as an estimate of what's achievable with sufficiently different algorithms.** It's a reasonable planning scenario for extensions of the existing architectures, but the evidence is insufficient to extrapolate to fundamentally new approaches.

There is also a fairly important opportunity in the chemistry resource formula that I think deserves more attention than Claude gives it: **the interaction between angle precision, factorization shape, and discontinuous qubit-width changes.**

## 1. What the new evidence actually establishes

The native chemistry verifier confirms the following breakdown:

| Component | Toffolis | Share of total |
|---|---:|---:|
| Inner state preparation | 17,333 | 56.8% |
| ANGLE lookups, both directions | 10,604 | 34.8% |
| Givens rotations | 2,250 | 7.4% |
| Outer state preparation | 317 | 1.0% |
| **Total** | **30,504** | **100%** |

The scored width is:

\[
Q=141+337=478.
\]

Those 337 temporary qubits are overwhelmingly associated with inner state preparation, not the ANGLE lookups. The current oracle architecture uses only 57 temporary qubits.

The verifier also confirms that your chemical reconstruction is extremely close to its acceptance boundary:

\[
\epsilon_2=19.9953045\text{ mHa},
\]

against a maximum of 20 mHa.

That leaves just **0.0047 mHa** of two-body error headroom. norm_cycle2_signed_verify

This establishes two things.

First, optimizing ANGLE workspace by itself is essentially pointless unless it changes the overall preparation architecture. The 337-qubit inner-preparation peak dominates.

Second, even small changes to the current quantized factorization can invalidate it. New structures need to be optimized for physical accuracy, not simply compressed from the existing one.

### What it does not establish

It does not establish that chemistry cannot improve 2×.

The metric charges

\[
C_{\mathrm{chem}}=\lambda_{\mathrm{eff}}TQ.
\]

It does not require using the current inner-preparation dimensions, the current angle encoding strategy, or a particular factorization shape.

Claude's projected maximum improvement comes from specific resource-model evaluations such as:

| Shape | Predicted chemistry improvement |
|---|---:|
| One block, 511 vectors, 350 columns | ~1.37× |
| Two blocks, 351 vectors and 175 columns each | ~1.45× |
| Same two-block shape with seven-bit angles | ~1.56× |

These estimates assume particular lookup costs and approximately unchanged normalization. None of the proposed shapes has been shown to produce an accepted quantized submission. claude_thoughts

**These are candidate architectures, not an exhaustive optimization of the cost model.**

That distinction is central to my disagreement with Claude's conclusion.

## 2. A specific chemistry counterexample to the proposed 1.6× ceiling

The round-4 report reveals that the host's inner-preparation circuit uses a discrete QROAM block-size parameter \(k_2\).

The relevant constraints are approximately:

\[
b_2=\lceil\log_2 B\rceil+11
\]

and

\[
(2^{k_2}-1)b_2\leq75b_{\mathrm{rot}}.
\]

Here \(B\) is the number of dictionary vectors per ensemble and \(b_{\mathrm{rot}}\) is the angle precision.

The host chooses \(k_2\) automatically. You cannot directly override it in a submission.

The consequence is that **small changes to angle precision or dictionary dimensions can change the number of temporary qubits by hundreds**. This is documented in the round-4 report, including an R24/B64/C32 example whose inner temporary width changes from 544 to 272 qubits between eight and seven bits. REPORT

### Consider R4/B180/C88

This is a four-block factorization with 180 vectors and 88 columns per block.

It has:

- 720 total dictionary vectors, almost identical to your existing 719.
- 352 total columns, almost identical to your existing 350.
- Approximately the same number of ANGLE lookup addresses.

Using the reported inner-preparation formulas:

| Resource | Eight-bit angles | Seven-bit angles |
|---|---:|---:|
| Inner word width \(b_2\) | 19 | 19 |
| Allowed \(k_2\) | 5 | 4 |
| QROAM block size | 32 | 16 |
| Inner temporary qubits | **607** | **304** |
| Inner-preparation Toffolis | ~3,412 | ~5,220 |

The seven-bit circuit uses more inner-preparation gates but approximately **303 fewer temporary qubits**.

That's a dramatic architectural difference.

If I retain roughly the incumbent ANGLE cost, use the reported Givens cost scaling, and estimate the remaining persistent storage from the incumbent, the seven-bit shape has approximately:

\[
T\approx18{,}100,\qquad Q\approx436.
\]

With unchanged effective normalization, that corresponds to a chemistry improvement of approximately:

\[
\boxed{1.85\times}
\]

If the smaller angle words and/or improved lookup construction reduce the ANGLE cost by 10–15%, the resource estimate improves to approximately **1.96–2.03×**.

**Important:** These are my conditional resource estimates, not native host measurements. The width estimate needs confirmation from the complete host model, the ANGLE cost is not measured, and no valid quantized factorization of this shape has been constructed.

In particular, seven-bit quantization may be extremely difficult at this shape.

Nevertheless, this illustrates the problem with Claude's conclusion.

The resource model does not obviously prevent 2× chemistry improvements. Claude evaluated a set of shapes that looked relatively achievable under previous numerical experiments, then effectively treated their best projected performance as an architectural ceiling.

That's not justified.

### An even stranger consequence: Adding vectors can reduce scored qubits

At eight-bit precision, increasing \(B\) from 256 to 257 increases \(b_2\) from 19 to 20.

This can force \(k_2\) from 5 down to 4, reducing the dominant difference-bank workspace term from

\[
31(19)=589
\]

to

\[
15(20)=300.
\]

So **adding one dictionary vector can remove approximately 289 temporary qubits from this component**.

The circuit would use additional Toffolis, but the reduction in \(Q\) could more than compensate.

This means the chemistry search must explicitly explore the discontinuities in the host resource model, including potentially padding a smaller factorization to a more favorable dimension.

It is particularly important not to assume that fewer vectors always produce a lower score.

### Does this restore my enthusiasm for joint diagonalization?

Partially.

The original attraction was that several small dictionaries might represent the same collection of matrices with considerably cheaper inner preparation.

The new information identifies the problem: a small dictionary can push the host into a larger QROAM block and drastically increase the width.

For example, my previous R10/B76 proposal is unattractive at eight-bit precision for precisely this reason.

**I would change the experiment from unrestricted joint diagonalization to resource-constrained multi-block factorization.**

Instead of asking which block structures reconstruct the tensor most accurately, first identify attractive regions of the actual resource landscape. Then use joint diagonalization, shared-basis clustering, and overcomplete dictionaries to try to fit those shapes.

R4/B180/C88 at seven bits is one candidate worth including in that resource-only search. It is not yet evidence that the corresponding physical fit exists.

## 3. DQI: The detailed report strengthens the case for combining the round-3 and round-4 architectures

This is probably the most important practical change to my recommendations.

Your round-3 and round-4 decoders have complementary strengths:

| Property | Round 3 | Round 4 |
|---|---:|---:|
| Architecture | Berlekamp–Massey | Compact Euclidean |
| Qubits | 1,429 | 961 |
| Toffolis | 645,755 | 1,019,866 |
| Phase errors | 8 | 380 |
| Repetition multiplier | 1.223 | 1.697 |
| DQI cost | 1.613T | 1.598T |

Round 4 saved 468 qubits but increased the gate count by approximately 58% and increased the repetition multiplier substantially. REPORT REPORT

### A particularly interesting target

Suppose you could achieve the following combination:

- Round-4 width: 961 qubits.
- Round-3 gate count: 645,755 Toffolis.
- Round-3 repetition multiplier: 1.223.

The resulting DQI cost would be approximately:

\[
C_{\mathrm{DQI}}=0.729\text{ trillion}.
\]

That's a **2.19× improvement over round 4**, without introducing a fundamentally cheaper decoder or phase-extraction algorithm.

It would raise the overall score from 9.37 to approximately 12.17 with the other tasks unchanged.

Of course, those parameters aren't independently selectable. The round-4 width reduction required changes that made arithmetic more expensive.

Nevertheless, this identifies a concrete architectural research problem:

> Can we preserve the compact storage layout of the round-4 Euclidean decoder while recovering the cheap arithmetic and low failure rate of the round-3 Berlekamp–Massey decoder?

I would give this substantial attention.

### Claude's exception-repair idea is good, but its explanation is not yet demonstrated

Claude attributes all 380 phase errors to the degree-one-quotient assumption.

That is plausible, given the history of the decoder, but the report does not contain a failure census establishing the causes of those 380 errors.

There's an intriguing coincidence:

- Your round-4 decoder fails approximately **0.371%** of shots.
- Gu et al.'s published half-degree Euclidean algorithm encounters an exceptional zero pivot in approximately **0.36%** of shots.

Gu et al. handle the first such exception and allow the much rarer second exception to remain uncorrected. [GitHub](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md?utm_source=chatgpt.com)

However, those are **different decoder algorithms and different Euclidean stages**. Similar failure frequencies do not establish that the same repair will work.

I would still make this the first DQI experiment, but its initial deliverable should be a failure classification, not a full implementation of Gu's correction mechanism.

Identify precisely which mathematical invariant first becomes invalid on failing inputs. Then determine whether a single exceptional-state record fixes most failures.

If it does, the opportunity is excellent.

Reducing the repetition multiplier from 1.697 to 1.223 gives approximately 1.39× before accounting for implementation overhead.

One thing I'd change in Claude's brief: its proposed limit of 40 additional qubits is unnecessarily restrictive as a hard rejection criterion.

At 1.05 million Toffolis and \(R=1.223\), the DQI decoder could grow from 961 to approximately **1,116 qubits** and still break even.

A compact correction remains preferable, but I wouldn't reject an effective repair merely because it needs 50–80 extra qubits. It could still be an important stepping stone toward a better hybrid decoder.

---

## 4. DQI: Claude's in-place transformation is worth investigating

Claude proposes replacing the evaluator polynomial \(\omega\) with

\[
F=\omega\Psi\bmod\sigma,
\]

where

\[
\Psi\sigma'\equiv1\pmod{\sigma}.
\]

This is a particularly good idea because, at each root \(\beta\) of \(\sigma\),

\[
F(\beta)=\frac{\omega(\beta)}{\sigma'(\beta)}.
\]

Thus, the error magnitude is available through evaluating \(F\), without performing a field division during the position scan.

The transformation is invertible whenever \(\Psi\) is invertible modulo \(\sigma\), as it is for a squarefree locator.

Consequently, an in-place reversible implementation is mathematically possible.

The question is whether it can be implemented **economically**.

A dense polynomial transformation with coefficients depending on the quantum locator is not a free change of basis. Performing the multiplication, preserving enough information to reverse it, and recovering the original representation could consume much of the proposed savings.

I would regard this as a promising research hypothesis rather than an established path to an 11-Toffoli position kernel.

There's an additional complication: the published Gu decoder uses a fixed classical mask, while the Kaggle verifier provides arbitrary classical measurement masks.

The published construction eliminates 71 scan positions using a correction specialized to the fixed mask. That simplification cannot automatically be reused without accounting for how the variable masks are processed. [Andi Gu](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html?utm_source=chatgpt.com) [Andi Gu](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html?utm_source=chatgpt.com)

**My recommendation:** Keep this as a separate high-risk DQI experiment, with a strict requirement to account for the complete reversible transformation, including its inverse and temporary storage.

I would investigate it before attempting an entirely new half-GCD decoder.

---

## 5. Chemistry: I agree with Claude about the quantizer, but not about its supposed precision floor

The chemistry experiments show a consistent and important problem:

| Experiment | Continuous fit | Quantized fit |
|---|---:|---:|
| Two-block B479 | ~13.0 mHa | ~25.0 mHa |
| Two-block B447 | ~14.8 mHa | ~27.8 mHa |
| Grouped R28/B64 | 17.02 mHa | 391.18 mHa at seven bits |

The control experiment is especially revealing: Claude's local quantizer failed to preserve a previously valid solution, producing approximately 22.7 mHa where earlier optimization had achieved around 19.6 mHa. REPORT

That is strong evidence that **the quantization pipeline itself is a significant limitation**.

But Claude then describes the approximately 8–9 mHa difference between continuous and quantized solutions as though it were an inherent eight-bit precision cost.

I disagree with that interpretation.

The best quantized representation need not be the nearest quantized approximation to the best continuous representation.

The optimizer might find a quite different discrete dictionary that exploits redundancy among vectors and compensation through the factorization weights.

In other words,

\[
\min_{\text{quantized }U,W}
\mathcal L(U,W)
\]

need not be well approximated by

\[
\mathcal L\left(
\operatorname{round}(U_{\text{float}}),
W_{\text{refit}}
\right).
\]

Even sophisticated local integer descent can remain trapped near the rounded continuous solution.

That leaves room for genuinely different approaches to quantized factorization, rather than simply improving the existing rounding procedure.

### What I would change

I would pursue two independent chemistry experiments.

**Experiment A: Repair and improve the current quantizer.**

This is Claude's suggestion, and I agree with it.

First reproduce the accepted incumbent. Then compare progressive freezing, rounding-error feedback, tree changes, and coupled integer-angle moves against the same floating-point checkpoints.

**Experiment B: Optimize a discrete factorization directly.**

Instead of starting with an accurate continuous dictionary and rounding it, construct candidate dictionaries directly from vectors representable by the permitted integer Givens angles.

Fit the weights against those vectors and iteratively modify the dictionary to improve the physical residual.

The two approaches explore different optimization landscapes.

I'd especially investigate the second approach for the seven-bit, four-block architectures, where reaching an entirely different quantized solution might be necessary.

### One other qualification

Claude's successful B479 continuous fit should not be confused with demonstrating that R2/B351/C175 is physically feasible.

The reported B479 refers to the number of vectors **per block**. The proposed R2/B351 structure is substantially smaller.

The direct R2/B351/C175 fit achieved approximately 21.1 mHa, which is already above the acceptance threshold *before quantization*. A larger-column R2/B351/C351 fit reached 17.2 mHa, but its resource cost is a different proposition. REPORT

Thus, even the proposed 1.45× chemistry improvement remains unproven on physical accuracy.

That is a legitimate reason for caution. It is not a reason to declare the architecture impossible.

---

## 6. I wouldn't abandon lookup optimization

Claude deprioritizes content-aware ANGLE optimization because it accounts for only 35% of the current chemistry Toffolis.

I think that is too dismissive.

At the current shape, reducing ANGLE Toffolis by 50% would lower total Toffolis by approximately 17.4%.

That alone is useful but not transformative.

However, suppose a smaller factorization reduces total Toffolis to approximately 21,840, including 10,380 ANGLE Toffolis, as Claude's two-block projection suggests.

Then reducing ANGLE Toffolis by 60% would lower the total to approximately 15,612.

At approximately 461 qubits and unchanged normalization, that combination would achieve roughly **2× chemistry improvement**.

Neither individual change does so.

This is why I would not treat the 35% share as a reason to stop investigating lookup algorithms.

The real question is whether the angle tables have exploitable structure.

Gu et al. demonstrated substantial lookup compression on a different chemistry benchmark using address scheduling and repeated data patterns. Their result does not guarantee that your FeMoCo tables are similarly compressible, but it provides a concrete precedent for content-aware circuit synthesis. [GitHub](https://github.com/andigu/uncompute/blob/main/README.md?utm_source=chatgpt.com)

A good experiment would measure repeated subtrees, bit dependencies, and Boolean-circuit complexity of your current accepted angle tables.

If the data appear nearly incompressible under those representations, terminate the experiment. Otherwise, investigate whether the numerical optimizer can encourage the kinds of patterns that the circuit compiler exploits.

---

# 7. How this changes the path to 17.71

The new information makes the trade-offs more concrete.

Using the actual round-4 local score of 9.3677, the following combinations would be required to match 17.71:

| Additional chemistry gain | Additional inverse gain | Required DQI gain | Target DQI Toffolis* |
|---|---:|---:|---:|
| 1.45× | 1.20× | 3.88× | ~364k |
| 1.60× | 1.30× | 3.25× | ~436k |
| 1.85× | 1.30× | 2.81× | ~504k |
| 2.00× | 1.30× | 2.60× | ~544k |

*Assumes 961 DQI qubits and a repetition multiplier of 1.223. These are derived targets, not predicted outcomes.*

This table is why I consider the chemistry investigation important.

If you achieve only 1.45× chemistry and 1.2× inverse, then reaching 17.71 demands an exceptionally good DQI decoder.

If chemistry can improve approximately 2×, the required DQI improvement becomes substantially more manageable.

And there's real research precedent for cheaper DQI arithmetic: the public Gu et al. construction uses approximately 440,000 Toffolis, though at a much larger width and with differences in the verification interface. [GitHub](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md?utm_source=chatgpt.com)

The goal is to combine low arithmetic cost with the compact storage already demonstrated by Codex.

---

# 8. Revised round-5 priorities

I would now organize the experiments as follows.

| Priority | Experiment | What would constitute an important result |
|---|---|---|
| **1** | DQI exceptional-case classification and repair | Reduce failures from 380 to approximately 8 without excessive additional width or gates |
| **2** | Chemistry exact resource-frontier search | Identify the best shapes across 7–9 bits, including QROAM width cliffs and actual ANGLE costs |
| **3** | Chemistry quantizer calibration and direct discrete fitting | Reproduce the incumbent, then materially reduce the gap between continuous and quantized physical error |
| **4** | DQI hybrid decoder | Combine round-3-style cheap arithmetic with round-4-style compact storage |
| **5** | Chemistry multi-block factorizations | Test resource-favorable shapes, including R4/B180/C88 at seven bits |
| **6** | DQI in-place \(F=\omega\Psi\bmod\sigma\) transformation | Establish a complete reversible implementation with a favorable total cost |
| **7** | Content-aware ANGLE synthesis | Measure whether existing and redesigned dictionaries admit substantially cheaper lookups |

I would keep inverse-square-root optimization running only if you have sufficient independent resources. It offers potential gains, but the first six experiments address larger opportunities.

**A critical change from my earlier recommendation:** I would not launch the multi-block fitting experiment before the exact resource-frontier scan.

The report explicitly acknowledges that no complete Pareto frontier using the current signed-gauge oracle architecture and host QROAM choices has been constructed. REPORT

That should be fixed first. Otherwise, an agent could spend many hours fitting a physically difficult shape that was never especially attractive under the actual scoring model.

For the numerical experiments, I would also avoid hard time-based or percentage-based rejection criteria until the important architecture-specific resource trade-offs are understood.

A failed optimizer, an unsuccessful quantizer, and a genuinely unfavorable cost bound are different outcomes.

---

# 9. Was Claude being defeatist?

I would separate its engineering judgment from its predictions.

### Where Claude is right

The current submission is highly optimized relative to its host reference.

Chemistry has an extremely tight physical-error margin. Its dominant workspace is determined by a discrete resource model, and the currently tested smaller factorizations have not survived quantization.

DQI is close to the mandatory syndrome-input width, so simply removing workspace will not produce another dramatic improvement.

Those are real obstacles.

I also think Claude's **12–14 range is reasonable as a planning scenario for incremental and moderately architectural improvements during one optimization round**.

### Where Claude overreaches

Its statement that chemistry tops out near 1.6× is not justified by the evidence presented.

Neither is an implied overall ceiling near 16.7.

The report describes incomplete resource searches, numerical optimizer failures, and quantization difficulties—not comprehensive lower bounds on the acceptable circuit families.

Moreover, round 4 already contains examples of Claude treating limitations of its own implementation as broader constraints, which were subsequently contradicted by Codex's verified work. REPORT

I wouldn't discard Claude's analysis because of those mistakes, but I would require unusually strong evidence before accepting another broad claim of infeasibility.

### My updated assessment

I am **less confident that the particular 2× chemistry architecture I originally suggested will work**, but more convinced that the chemistry cost landscape has not been adequately explored.

The seven-bit width discontinuity is especially interesting. It creates configurations with attractive resource estimates that have not been subjected to serious quantization-aware optimization.

Meanwhile, DQI contains two complementary partially successful architectures and at least one promising unfinished phase-extraction approach.

That is enough to justify targeted research.

I would distinguish three statements:

- **A score of 12–14 is a plausible result from extending existing methods.**
- **A score above 17 requires substantial architectural progress and is genuinely difficult.**
- **Nothing in these documents demonstrates that 17+ is beyond reach.**

I agree with the first two. I see no convincing evidence for the third statement's opposite.

**For round 5, my strongest recommendation is to pursue DQI correctness repair and the exact chemistry resource-frontier search first. Then use those results to select the more radical implementations.**

In particular, I would make Codex investigate the **seven-bit, multi-block chemistry width cliffs**, rather than assuming the best available chemistry architecture must look like the existing 719-vector dictionary or Claude's two-block proposal.

That's a concrete, previously underexplored direction where the reported cost model suggests substantial potential—and where a relatively inexpensive resource-only experiment can establish whether deeper numerical research is worthwhile.