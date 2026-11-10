# Research assessment: Breaking past a score of 9.38

## 1. The central problem is architectural, not incremental

After reviewing the eight supplied documents and the published research, I agree that **round 5 should prioritize a small number of architectural experiments rather than another open-ended optimization run**.

Your existing agents have already made substantial progress:

| Task | Original cost | Round 4 cost | Improvement |
|---|---:|---:|---:|
| Inverse square root | 2.077M | 182,957 | 11.35× |
| DQI decoding | 24.439T | 1.598T | 15.29× |
| Quantum chemistry | 3.356B | 708.55M | 4.74× |

The important observation is that your progress has come in discontinuous jumps.

- Round 1's major breakthrough was replacing the inverse-square-root algorithm with direct piecewise polynomial approximation.
- Round 2's was replacing the original synchronized Euclidean decoder.
- Round 3's was replacing the DQI implementation with a syndrome-consuming Berlekamp–Massey architecture and introducing a shared chemistry dictionary.
- Round 4 mostly improved those existing architectures rather than replacing them.

That history supports your concern about local optima, although it doesn't establish that the current architectures are near their global limits. optimization_iteration_1 optimization_iteration_2 optimization_iteration_3 optimization_iteration_4

### How big a breakthrough is actually needed?

The competition uses the geometric mean of the three task improvements. Consequently,

\[
S_{\text{new}}=9.38
\left(
F_{\text{inverse}}F_{\text{DQI}}F_{\text{chemistry}}
\right)^{1/3},
\]

where each \(F\) represents an additional reduction in that task's cost.

To reach **17.0**, you need:

\[
F_{\text{inverse}}F_{\text{DQI}}F_{\text{chemistry}}
\geq 5.95.
\]

To reach the current leader's reported **17.71**, you need approximately **6.73×**.

There are several ways to get there:

| Additional improvements | Projected score |
|---|---:|
| 6× DQI, other tasks unchanged | 17.05 |
| 3× DQI, 2× chemistry | 17.05 |
| 2× DQI, 2× chemistry, 1.5× inverse | 17.05 |
| 2× improvement in all three | 18.76 |

This changes how I would prioritize research.

A 10% improvement to one task increases the overall score by only about 3.2%. Even a 50% improvement to one task yields only about 14% overall.

**The target is attainable mathematically without a 6× breakthrough in any single task**, but at least one substantially new algorithmic approach probably needs to work.

### Where the mathematical leverage lies

DQI is particularly attractive because its cost scales as

\[
C_{\mathrm{DQI}}=R(\epsilon_+)\,Q^2T.
\]

However, your round-4 implementation already uses 961 qubits, of which **840 are mandatory syndrome-input qubits**. Further width reduction has limited room.

For example, even reducing DQI width from 961 to the theoretical input-only floor of 840 would improve the task cost by just 1.31× if everything else stayed constant.

To obtain a 6× DQI-only improvement, an implementation using 900 qubits and the near-ideal repetition multiplier of 1.19 would need approximately **276,000 Toffolis**, compared with your current 1.02 million.

That is an extremely aggressive target, especially because the public Berlekamp–Massey implementation already requires approximately 440,000 Toffolis at a much larger width.   [GitHub](https://github.com/andigu/uncompute/blob/main/README.md?utm_source=chatgpt.com)

I would therefore pursue **DQI architecture changes and chemistry architecture changes simultaneously**, with inverse-square-root optimization as a secondary opportunity.

## 2. DQI: A cheaper phase-extraction architecture

**Priority: High | Potential: Moderate individually, substantial when combined with decoder changes**

One useful finding from the public research is that the phase-extraction kernel can be considerably cheaper than the one currently used in your submission.

The published `andigu/uncompute` Reed–Solomon decoder uses a **23-Toffoli kernel per field position**, consisting of:

- 11 Toffolis to test whether the position is a root.
- 12 Toffolis to compute the required trace pairing and apply the phase.

Your round-4 implementation reports an 88-Toffoli point kernel. It is therefore worth investigating whether the published construction can be adapted to your narrower decoder.  

**Research:** [Andi Gu's Reed–Solomon decoder, with annotated implementation](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html?utm_source=chatgpt.com)

### How the 23-Toffoli construction works

Ordinary Reed–Solomon decoding uses Forney's formula:

\[
e_j=\frac{\omega(\beta_j)}{\sigma'(\beta_j)}
\]

at roots \(\beta_j\) of the error-locator polynomial \(\sigma\).

A naive implementation computes an inverse at each position. That's expensive because there are 4,095 possible positions.

The published construction instead computes a polynomial \(\Psi\) satisfying

\[
\Psi(x)\sigma'(x)\equiv1\pmod{\sigma(x)}.
\]

Consequently, at each root,

\[
e_j=\omega(\beta_j)\Psi(\beta_j).
\]

It then uses two additional observations:

1. Squaring over a binary extension field is linear, allowing the inverse-derivative problem to be reduced to a roughly half-degree Euclidean problem.
2. The required phase depends on a trace pairing, not on the complete error magnitude. The trace pairing can be compiled directly as a bilinear Boolean form.

The resulting implementation avoids constructing a complete 12-bit error value at every position.

The published implementation also uses a linear correction derived from the syndromes to eliminate 71 positions from the expensive scan. [Andi Gu](https://andigu.github.io/projects/uncompute/dqi-reed-solomon.html?utm_source=chatgpt.com)

### Why this matters for round 5

Replacing an 88-Toffoli kernel with a 23-Toffoli kernel across 4,095 positions would nominally save:

\[
4095(88-23)=266{,}175\text{ Toffolis}.
\]

Against your current 1.02 million Toffolis, that's approximately a 26% reduction in the gate count of the complete DQI circuit.

**However, this is not a free 266,000-Toffoli improvement.**

The cheaper kernel needs a suitable representation of the locator, evaluator, and inverse-derivative polynomial. Constructing those objects and subsequently uncomputing them has a cost. The published decoder also uses substantially more qubits than your current implementation.

A simple transplant could therefore score worse.

The research opportunity is to combine:

- Your 961-qubit storage architecture.
- The published trace-pairing phase extraction.
- A cheaper way to obtain the inverse-derivative information.
- Better handling of exceptional decoding paths to reduce the repetition multiplier.

This is sufficiently concrete that I would give Codex a dedicated experiment rather than asking it to optimize DQI generally.

**Initial experiment:** Make the existing decoder expose whatever polynomial data the new kernel needs, implement only the forward phase path, and calculate the complete projected cost of generating and removing that data. Do not rebuild the entire decoder until the projection beats the incumbent.

## 3. Chemistry: Replace the shared dictionary with jointly diagonalized blocks

**Priority: Very high | Risk: High | Potential: Large**

This is the most interesting chemistry direction I found. It is also a way to revisit the failed multi-block experiments without repeating the same optimization approach.

Your current factorization uses a large shared dictionary, with 719 vectors and 350 columns. The architecture has excellent expressiveness, but the large dictionary drives state-preparation and lookup resources.

One alternative is to use **multiple small dictionaries instead of one enormous dictionary**.

The important distinction is *how those dictionaries are constructed*.

### A constructive way to generate a multi-block factorization

The chemistry Hamiltonian's two-electron tensor can be represented as

\[
H_2\approx\sum_{c=1}^{C}X_c\otimes X_c,
\]

where each \(X_c\) is a symmetric \(76\times76\) matrix.

Your current dictionary expresses these matrices as

\[
X_c=\sum_{b=1}^{719}w_{bc}u_bu_b^\mathsf T.
\]

Now consider an alternative construction.

Because \(X_c\) is symmetric, it has an eigendecomposition:

\[
X_c=V_cD_cV_c^\mathsf T.
\]

That means each individual \(X_c\) can be represented using **at most 76 orthogonal vectors**, rather than 719 shared vectors.

This immediately gives an exact reparameterization of the *current represented two-body tensor*, before angle quantization: assign each \(X_c\) its own dictionary of eigenvectors.

Unfortunately, doing that independently for all 350 matrices would require up to \(350\times76\) dictionary vectors. The lookup cost would become enormous.

But what if groups of \(X_c\) matrices share approximately the same eigenvectors?

Then a group could use one shared orthogonal basis:

\[
X_c\approx V_\rho D_c V_\rho^\mathsf T
\qquad(c\in G_\rho).
\]

The resulting factorization has the same mathematical form as the competition's allowed multi-block factorization.

For example, a candidate with

\[
(R,B,C)=(10,76,35)
\]

would contain 760 dictionary vectors and 350 coefficient columns across its blocks.

The total lookup-table size would be comparable to your existing 719-vector dictionary, but **the inner state preparation would operate on blocks of 76 vectors instead of 719**.

Whether this actually reduces total scored width depends on the verifier's resource formula. That is measurable before building the full submission.

### The new algorithm: Joint approximate diagonalization

Finding a basis that approximately diagonalizes several matrices is a well-established numerical optimization problem.

Given a group of symmetric matrices, find an orthogonal \(V\) minimizing

\[
\mathcal L(V)=
\sum_{c\in G}
\left\|
\operatorname{offdiag}(V^\mathsf T X_cV)
\right\|_F^2.
\]

This is called **joint approximate diagonalization**.

Relevant research:

- [Randomized Joint Diagonalization of Symmetric Matrices](https://epubs.siam.org/doi/pdf/10.1137/22M1541265?download=true&utm_source=chatgpt.com) — He and Kressner, SIAM.
- [Jacobi Angles for Simultaneous Diagonalization](https://epubs.siam.org/doi/abs/10.1137/S0895479893259546?utm_source=chatgpt.com) — Cardoso and Souloumiac.
- [pyRiemann's existing JADE implementation](https://pyriemann.readthedocs.io/en/latest/generated/pyriemann.geometry.ajd.jade.html?utm_source=chatgpt.com).

The crucial caveat is that arbitrary symmetric matrices do not share an eigenbasis. If the 350 matrices are strongly noncommuting, a small number of blocks may be unable to approximate them within the competition's 20 mHa tolerance.

That is why the first experiment should be a numerical feasibility test, not circuit generation.

### A second degree of freedom makes this more interesting

There is another transformation the current work has already used, but primarily for improving normalization.

Suppose the collection of matrices is arranged into a vector:

\[
X=(X_1,\ldots,X_C).
\]

For any orthogonal matrix \(O\),

\[
\widetilde X_c=\sum_jX_jO_{jc}
\]

preserves the represented two-body tensor exactly:

\[
\sum_c\widetilde X_c\otimes\widetilde X_c
=
\sum_cX_c\otimes X_c.
\]

The identity offsets can be rotated correspondingly.

Thus, before trying to group and jointly diagonalize these matrices, we can **rotate the entire collection of generators to make it easier to group**.

This leads to a more interesting optimization:

1. Orthogonally mix the 350 generator matrices.
2. Cluster them into groups based on how well they can share an eigenbasis.
3. Jointly diagonalize each group.
4. Construct the corresponding multi-block factorization.
5. Refit coefficients and identity offsets.
6. Optimize the actual eight-bit angle representation and resulting resource cost.

The mode rotation is not guaranteed to improve joint diagonalizability, but it substantially enlarges the search space beyond simply partitioning the current dictionary vectors.

### Why the previous failures don't settle this

Round 3 produced a grouped floating-point fit with approximately **14.07 mHa** two-body error. It subsequently deteriorated to approximately 31 mHa after eight-bit quantization.

Round 4's Claude experiments also found warm-started smaller dictionaries that achieved approximately 13 mHa before quantization, but failed after quantization.

Those are important results: they demonstrate that at least some alternative factorization shapes have promising continuous solutions. They do not establish that the shapes are impossible at eight-bit precision. optimization_iteration_3 optimization_iteration_4

AGY's disjoint-block attempt stalled at 479.36 mHa. I would not interpret that as strong evidence against jointly diagonalized blocks.

A naive partition of the original dictionary generally destroys cross-block contributions to the tensor. Grouping the **generator matrices** and finding new bases is a different operation.

**My proposed round-5 experiment:** Start with the accepted 719-vector solution, construct its 350 generator matrices, and investigate the approximation error attainable with jointly diagonalized groups before attempting quantization.

Try different combinations of \(R\), \(B\), and \(C\), with an initial emphasis on \(B\leq76\).

If the continuous solutions are promising, invest in quantization-aware fitting. If their best achievable tensor residuals are already well above 20 mHa, terminate the approach early.

This is a genuine architectural search with the potential to change the qubit-width bottleneck.

## 4. DQI: Fix rare decoding failures instead of accepting them

**Priority: High | Risk: Low–Medium | Potential: Moderate**

This is less radical than changing the decoder completely, but the expected payoff is unusually attractive.

Your round-4 decoder reports:

| Metric | Value |
|---|---:|
| Qubits | 961 |
| Toffolis | ~1,020,000 |
| Phase errors | 380 / 102,400 |
| Repetition multiplier | 1.69685 |

The repetition multiplier is currently quite expensive.

A decoder with zero observed failures at the same default shot count receives a multiplier of approximately 1.19126 rather than 1.69685.

At unchanged width and gate count, eliminating the failures would reduce the DQI cost by approximately **30%**, yielding a 1.425× DQI improvement.

You could even increase Toffolis significantly and still come out ahead.

### Suggested approach: Exceptional-path reconstruction

The public `uncompute` project provides useful precedent. Its Berlekamp–Massey implementation reconstructs most overwritten information from the surviving registers, but reserves a small amount of storage for exceptional cases.

The researchers evaluated 403 million sampled decodings and found that a single 19-qubit exceptional-history entry was sufficient almost all the time. Reserving space for the theoretical maximum would have required 665 qubits. [Anthropic](https://www.anthropic.com/research/designing-efficient-quantum-circuits?trk=public_post_comment-text\&utm_source=chatgpt.com)

Your exceptional conditions are not necessarily identical. However, the underlying strategy transfers:

1. Instrument the current decoder to identify the exact causes of the 380 phase failures.
2. Measure how frequently each exceptional condition occurs.
3. Determine the minimum additional reversible information needed to handle each condition.
4. Add targeted correction machinery rather than making the entire decoder general-purpose.
5. Re-evaluate the complete \(R(\epsilon_+)Q^2T\) objective.

One caution: **a rare quantum-controlled branch does not automatically have a cheap average Toffoli count**. The verifier counts executed quantum gates even if their quantum controls happen to be zero. A fallback must therefore be inexpensive as a circuit, not merely infrequently activated.

I would test this independently of the 23-Toffoli phase kernel.

The two improvements are complementary.

## 5. DQI: Replace sequential decoding with recursive rational reconstruction

**Priority: Medium | Risk: Very high | Potential: Large**

This is a more speculative architecture change.

Both your existing implementation and the public Berlekamp–Massey decoder perform sequences of polynomial updates.

That creates two significant costs: the arithmetic needed to find the error locator and the information that must survive so those updates can be reversed.

An alternative is to treat decoding as a polynomial rational-reconstruction problem and solve it recursively.

The mathematical relationship is the Reed–Solomon key equation:

\[
\sigma(x)S(x)\equiv\omega(x)\pmod{x^{70}}.
\]

Instead of 70 sequential Berlekamp–Massey updates, one could investigate divide-and-conquer Padé approximation, half-GCD, or minimal approximant-basis algorithms.

Recent work gives some reason to investigate this:

- [Refined complexity bounds for rational reconstruction and XGCD through Padé approximants and Cauchy interpolants](https://arxiv.org/abs/2609.21051?utm_source=chatgpt.com) (September 2026).
- [Optimizing the half-GCD algorithm](https://researchportal.ip-paris.fr/en/publications/optimizing-the-half-gcd-algorithm/?utm_source=chatgpt.com) (2026).
- [FLINT's Reed–Solomon decoding implementation](https://fredrikj.net/fdoc/html/tutorials/algebra/reed-solomon.html?utm_source=chatgpt.com), which uses Berlekamp–Massey and faster polynomial algorithms.

These approaches reduce asymptotic arithmetic complexity, but that does **not** automatically imply a smaller reversible circuit at degree 70.

Their intermediate polynomial products and recursion histories may outweigh the savings.

### Why this is worth investigating anyway

Round 4 already produced a component-level result involving Wronskian padding: adapter work fell from 28,069 to 7,140 Toffolis per direction, with component tests passing.

A complete compatible decoder was never constructed. optimization_iteration_4

That makes this an unfinished lead, not a failed algorithm.

My suggested experiment is to give Codex a separate branch specifically for the Wronskian/rational-reconstruction architecture.

The deliverables should be:

- A classical reference implementation reproducing the required locator/evaluator.
- A reversible polynomial-arithmetic implementation with measured forward and reverse resources.
- A complete qubit-lifetime schedule.
- A projected end-to-end DQI cost, including the phase kernel and repetition penalty.

**Gate-budget target:** With 1,000 qubits and a repetition multiplier of 1.20, approximately 444,000 total Toffolis would deliver a 3× DQI improvement over round 4.

That is ambitious, but it gives the agent a concrete architectural target rather than an arbitrary leaderboard goal.

## 6. Chemistry: Quantization-aware fitting needs to become the primary optimizer

**Priority: Very high | Risk: Medium | Potential: Large**

There is a recurring pattern in the chemistry experiments that I think deserves more attention.

Your agents have repeatedly produced promising continuous fits, followed by substantial deterioration when converting vectors to Givens angle words.

Examples include:

| Experiment | Before quantization | After quantization |
|---|---:|---:|
| Round 3 grouped fit | 14.07 mHa | ~31 mHa |
| Round 4 Claude, B479 | ~13 mHa | ~25 mHa |
| Round 4 Claude, B447 | Not separately established | ~27.8 mHa |
| Round 4 AGY, 28-group fit | 17.02 mHa | 391.18 mHa at seven bits |

These are not necessarily failures of the underlying factorization architectures. They demonstrate failures of particular combinations of factorization, angle representation, and quantization procedure. optimization_iteration_3 optimization_iteration_4

### Proposed change

Instead of treating quantization as a conversion applied after fitting, optimize the actual discrete variables from the beginning of the final optimization stage.

The optimization variables should include:

\[
\{\,\text{integer angle words},\,W,\,H2_I,\,\text{tree structure}\,\}.
\]

The objective should evaluate the physical vectors decoded from those exact words, not nearby floating-point vectors.

One possible algorithm is:

1. Start from a good continuous factorization.
2. Quantize its angles.
3. Refit all linear weights and identity offsets against those fixed, quantized vectors.
4. Estimate which individual angle-word changes would improve the true reconstruction residual.
5. Apply a small collection of promising integer changes.
6. Refit the weights and offsets again.
7. Periodically propose coupled changes to multiple angles, avoiding a search confined to single-coordinate local optima.

I would also test annealed multicoordinate search and small-block discrete optimization, rather than relying entirely on straight-through estimators. Round 3 already reports that straight-through angle fitting failed to produce a better selected candidate.

The distinguishing feature would be **joint discrete moves followed by an accurate weight refit**, not simply another version of conventional quantization-aware training.

AGY's reported analytical offset solver, reducing an approximately 11-minute optimization to 4.4 seconds, could be especially useful here. That result should be reproduced and benchmarked before relying on it. optimization_iteration_4

### Optimize the tree too

The Givens tree is a particularly interesting variable because it changes how the same physical vector is represented by a sequence of quantized angles.

The competition explicitly allows changing the tree, and its 75 edges still represent 76-dimensional unit vectors.

A different tree can change:

- How angle quantization errors propagate into physical vector components.
- The prevalence of near-zero and otherwise simple angle words.
- Which angle bits are correlated or redundant across dictionary entries.
- The compressibility of the resulting lookup tables.

I would compare a standard balanced tree against trees constructed from amplitude statistics, orbital correlations, and direct optimization of quantized reconstruction error.

The result would then be evaluated with the actual oracle builder, because a better numerical tree need not have cheaper lookup circuits.  

### Why this might be disproportionately valuable

Suppose a multi-block architecture achieves a 478-to-260 qubit reduction while retaining approximately 28,000 Toffolis and the current effective normalization.

Its chemistry cost would be roughly half the current value.

That alone, combined with a 3× DQI improvement, would put the combined score near 17.

Those numbers are **design targets, not predicted performance**. Nevertheless, they show why getting a genuinely different chemistry factorization through quantization is potentially more valuable than another 1–2% improvement to the existing dictionary.

## 7. Inverse square root: Revisit piecewise-linear approximation

**Priority: Medium | Risk: Medium | Potential: Moderate**

I found a useful precedent for a different inverse-square-root architecture.

The published `uncompute` Coulomb-potential implementation replaces piecewise quadratics with a larger table of straight lines.

Its design uses:

- Adaptive interpolation regions.
- A short slope rather than a full quadratic coefficient.
- Fixed-point fitting directly against acceptable output intervals.
- Measurement-based deletion of lookup data.
- A specialized phase correction that costs 17 Toffolis instead of replaying the complete 126-Toffoli lookup.

This is a different benchmark, with a different accuracy contract and input construction, so its resource counts do not transfer directly to Kaggle's inverse-square-root task. [Andi Gu](https://andigu.github.io/projects/uncompute/coulomb-potential.html?utm_source=chatgpt.com)

### Proposed experiment

Your existing design uses 384 normalized intervals and a quadratic approximation.

I would construct an independent family of piecewise-linear circuits using approximately 512–2,048 intervals, with nonuniform segment widths.

The idea is to exchange more lookup entries for significantly less arithmetic.

A good search procedure would optimize:

\[
(\text{segment boundaries},\text{slope precision},
\text{intercept precision},\text{workspace},\text{phase margin})
\]

jointly.

The most important detail is that fitting should use the **actual integer arithmetic and acceptance windows**, rather than first fitting ideal real-valued straight lines.

Then measure the complete circuit, including lookup deletion and rotation synthesis.

The target is to beat your current approximately 183,000 raw task cost, not simply its Toffoli count.

In principle, a narrower design with more Toffolis could still win because inverse-square-root cost is proportional to qubit count.

I would regard a 1.3–1.5× improvement as a worthwhile exploratory target, not a reliable expectation.

This is less likely to produce the largest single breakthrough, but it provides a plausible third contributor to a score above 17.

## 8. Chemistry: Compress lookup tables by optimizing the data representation

**Priority: Medium | Risk: Medium | Potential: Moderate**

One additional direction deserves a separate experiment: optimizing the mathematical representation specifically to make the QROM circuits cheaper.

Most of your chemistry fitting has focused on the physical tensor error and normalization. Round 4 made progress through signed-angle gauges, which rendered two payload bits dependent.

I suspect there is additional room for *representation-aware fitting*.

The idea is to treat the lookup circuit's Boolean complexity as part of the numerical optimization objective.

For example, consider two dictionaries that achieve essentially identical reconstruction error and normalization. One has angle words distributed almost randomly across addresses; the other has many correlated or repeated patterns.

The latter may admit a cheaper Boolean circuit even if both have identical dimensions and precision.

The public `uncompute` chemistry work provides a related example. Its Fe₂S₂ circuit substantially reduces lookup cost by reordering entries, exploiting repeated data, and pruning lookup-tree branches. That construction addresses a different chemistry problem, but demonstrates the value of optimizing data layout and Boolean-circuit structure together. [Andi Gu](https://andigu.github.io/projects/uncompute/chemistry-fe2s2.html?utm_source=chatgpt.com)

### Specific ideas

For the Kaggle ANGLE tables, I would investigate three transformations:

**First, dictionary permutation.** The dictionary vectors may be reordered as long as their coefficient rows are reordered consistently. Search for arrangements that minimize measured lookup-circuit cost.

**Second, angle gauge optimization.** Extend round 4's signed-angle transformations into a systematic search over equivalent angle encodings. Prefer encodings with correlated bit patterns and simpler Boolean functions.

**Third, content-aware lookup synthesis.** Instead of paying for a generic lookup tree regardless of the data, use a shared multi-output Boolean representation that exploits identical subtrees, repeated bit patterns, and constant outputs.

The third idea has the most implementation risk.

The key experiment is to compare measured Toffoli counts for the *existing accepted angle words*, so that physical accuracy and normalization remain exactly unchanged.

Only if that produces savings should Codex attempt joint optimization of dictionary coefficients and QROM complexity.

This would isolate algorithmic lookup improvements from numerical-factorization changes.

---

## 9. Approaches I would deprioritize

Your previous experiments and the verifier's cost formulas rule out several superficially promising directions as primary round-5 projects.

| Approach | Reason for deprioritizing |
|---|---|
| Further chemistry normalization-only optimization | Current effective normalization is 48.594; the reported lower bound is 46.786. That leaves only about 3.7% cost headroom from normalization alone. |
| Further substantial DQI width reduction | At 961 qubits, the current decoder is already close to the mandatory 840-qubit input width. |
| Another generic GF(4096) multiplier search | The existing 42-Toffoli multiplier has survived extensive searches. An improvement would help, but it would not address the largest architectural costs. |
| Simple chemistry dictionary truncation | Multiple attempts reduced dimensions while increasing normalization or violating the physical-error threshold. |
| Naively lowering chemistry angle precision | Seven-bit quantization has repeatedly produced unacceptable errors. A different architecture or discrete optimization is needed first. |
| Replacing the DQI position scan with a sparse root list | Extracting 35 roots is not enough. The decoder must still evaluate an arbitrary classical mask at those locations, which introduces a potentially expensive lookup problem. |

I would not permanently abandon these directions. They simply have weaker evidence for producing the kind of discontinuous improvement you need.

One especially important distinction is that the 23-Toffoli phase kernel is already quite economical for its particular pointwise design. Its zero test costs 11 gates, and the nondegenerate 12-bit trace pairing has bilinear rank 12.

That suggests concentrating subsequent DQI work on **producing the polynomial data cheaply**, rather than continually trying to shave another gate from the point kernel.

---

# 10. My recommended round-5 experiment portfolio

I would organize round 5 around six separate experiments, with each preserving the current 9.38 submission as a fallback.

| Priority | Experiment | Immediate objective | Principal uncertainty |
|---|---|---|---|
| **1** | Chemistry: jointly diagonalized blocks | Find a low-width multi-block shape with <20 mHa quantized error | Can the generator matrices share sufficiently compact bases? |
| **2** | DQI: 23-Toffoli phase kernel | Combine cheap trace pairing with the narrow decoder | Does the inverse-derivative preparation erase the savings? |
| **3** | Chemistry: quantization-aware refitting | Recover valid eight-bit versions of B479 or grouped fits | Can discrete optimization close the 5–15 mHa quantization gap? |
| **4** | DQI: exceptional-path repair | Reduce the repetition multiplier toward 1.19 | How expensive is correcting the actual failure conditions? |
| **5** | DQI: Wronskian/rational reconstruction | Complete the abandoned alternative decoder | Can its arithmetic and storage compete end-to-end? |
| **6** | Inverse: adaptive piecewise-linear model | Reduce arithmetic and live workspace simultaneously | Can more lookups offset the tighter approximation requirements? |

I would give the first three experiments the most attention.

The Wronskian architecture is potentially very interesting, but its current evidence is primarily component-level. I would not assign it unlimited development time unless its projected complete resources are promising.

## Exact performance targets

For DQI, these are useful targets for Codex.

Assuming a repetition multiplier of 1.20:

| DQI qubits | Toffolis needed for 2× improvement | Toffolis needed for 3× improvement |
|---|---:|---:|
| 900 | 822,000 | 548,000 |
| 961 | 721,000 | 481,000 |
| 1,000 | 666,000 | 444,000 |
| 1,100 | 550,000 | 367,000 |
| 1,200 | 462,000 | 308,000 |

These illustrate why simply adapting the public 440,000-Toffoli decoder is not enough if it retains its approximately 1,837-qubit width.

The major opportunity comes from combining its arithmetic with your substantially smaller storage architecture.

For chemistry, the equivalent target is simpler.

Your current cost is approximately:

\[
48.594\times30{,}504\times478
=708.55\text{ million}.
\]

Assuming effective normalization remains unchanged, a **2× chemistry improvement** requires:

\[
Q_{\mathrm{chem}}T_{\mathrm{chem}}
\lesssim7.29\text{ million}.
\]

For example, 260 qubits and 28,000 Toffolis would satisfy this target.

That gives the multi-block experiment a concrete numerical goal.

### Three possible routes to 17+

| Scenario | DQI gain | Chemistry gain | Inverse gain | Overall score |
|---|---:|---:|---:|---:|
| A: DQI-heavy | 3.0× | 2.0× | 1.0× | **17.04** |
| B: Balanced | 2.5× | 1.8× | 1.35× | **17.12** |
| C: DQI plus modest gains elsewhere | 2.8× | 1.7× | 1.3× | **17.22** |

All three are speculative combinations, not forecasts.

I find scenario B particularly interesting because it doesn't require any one task to achieve an extraordinary improvement.

---

## 11. How I would structure the agent's work differently

Your concern about agents becoming trapped in local optima suggests changing more than the prompts.

I would change **what constitutes a successful experimental run**.

Currently, a goal such as `/goal score 14.0` incentivizes the agent to search for improvements anywhere. As progress slows, it is natural for the agent to prefer modifications that have already worked.

For radical research, that is counterproductive.

I would instead define a successful experiment as one that either produces a competitive architecture **or convincingly eliminates a specific hypothesis**.

### Each experiment should have three stages

**Stage A: Mathematical and numerical feasibility**

Before circuit construction, demonstrate that the proposed algorithm does what it claims.

For chemistry, this means testing the factorization and the exact quantized physical residual.

For DQI, it means constructing an independent classical implementation, checking the required phase against known error vectors, and measuring exceptional-path frequencies.

**Stage B: Component circuit feasibility**

Build the expensive new component with the correct reversible interface.

Measure actual Toffolis, maximum qubit lifetimes, and uncomputation costs.

In particular, avoid crediting a proposed gain before the component can be cleaned up reversibly.

**Stage C: Complete submission feasibility**

Only after the new architecture demonstrates promising component costs should Codex integrate it into a complete candidate and run the unchanged native verifier.

This is where the actual repetition multiplier, persistent qubits, analytic chemistry costs, and synthesis margins get accounted for.

### Require failure diagnostics

I would explicitly prohibit reporting that an idea is mathematically infeasible merely because an optimizer failed.

For every failed experiment, the agent should distinguish:

- A proved limitation of the algorithm.
- A counterexample or measured lower bound for the particular implementation.
- Numerical optimization failing to converge.
- Quantization destroying a valid floating-point solution.
- A circuit implementation exceeding the resource budget.
- A candidate that passes verification but scores worse.

Those have very different implications for future research.

Your round-4 reports contain a particularly relevant example: Claude subsequently retracted several unsupported conclusions about resource limits and the overall score ceiling. Codex's improvements were direct counterexamples. optimization_iteration_4

### Add one control that I think is essential

Before attempting new chemistry architectures, require the agent to demonstrate that its current fitting, quantization, and scoring pipeline reproduces the accepted incumbent.

Round 4 reported that one quantizer degraded an inherited-basin control to 22.7 mHa, despite the previous solution achieving approximately 19.6 mHa.

That is a serious warning sign.

If the transformation pipeline cannot preserve an already accepted solution, negative results on new architectures are difficult to interpret.

I would make successful round-trip reproduction a prerequisite for further chemistry experiments.

### Preserve experimental independence

For the DQI decoder, use the same representative error vectors to compare alternative algorithms during development. Then use independent full-verification runs and salts for candidate selection.

For chemistry, compare alternatives against the identical checksum-validated reference tensors and exact quantization semantics.

For inverse square root, reduced-domain tests are useful for debugging, but the final 3,133,452-input verification is mandatory.

The official scoring rules make these distinctions important: DQI samples and statistical confidence affect cost, while the other two tasks have hard correctness and physical-accuracy requirements. verifiers

---

# 12. Final assessment

I do **not** think your current results establish that a score near 9.4 is a fundamental ceiling for ordinary development agents.

They establish something narrower: several agents have converged on a set of increasingly optimized architectures, and incremental modifications to those architectures are producing diminishing returns.

That is exactly when an architecture-focused research campaign becomes attractive.

Of the ideas above, my strongest recommendations are:

**1. Revisit multi-block chemistry through joint diagonalization and orthogonal generator mixing.** This is my highest-risk, potentially highest-reward suggestion. It attacks the structure responsible for the chemistry width bottleneck, not just the current numerical parameters.

**2. Integrate the 23-Toffoli DQI trace kernel with the narrow decoder, while separately investigating rare-error correction.** These are grounded in working constructions and offer a more predictable route to improvement.

**3. Treat chemistry quantization as part of the primary optimization problem.** The gap between successful floating-point fits and failed quantized fits is too consistent to ignore. In particular, I would not consider the B479 and grouped architectures thoroughly explored.

**4. Finish evaluating the Wronskian alternative.** There is positive component-level evidence, but no complete circuit. It deserves a targeted experiment rather than a speculative claim about its score.

There is also an important limitation to my assessment: the supplied documents are summaries. They include useful measurements and detailed descriptions, but not the complete round-4 agent workspaces, source implementations, or intermediate candidate artifacts. I can identify which conclusions are unsupported, but cannot independently audit every AGY or Claude implementation from these materials alone.

**My preferred strategy would be to make round 5 a contest between fundamentally different representations of the chemistry tensor and the DQI decoding problem.**

The strongest path to 17+ appears to be some combination of substantially cheaper DQI decoding and a chemistry factorization with materially smaller state-preparation width.

I would prioritize discovering those architectures over extracting another few percent from the current solutions.