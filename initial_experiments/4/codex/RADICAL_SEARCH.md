# Active goal and current research

The user explicitly reactivated the14.0 goal, asks for radical changes, reports
a16.6 competitor result, and says to continue until achieved or interrupted.
The goal API now reports ACTIVE. Do not pause or mark complete below14.0.
No subagents have been authorized. All new artifacts remain in4/codex.

The current audited archive scores **9.332656488711933**, ZIP SHA256
81c3ea50eb2a3e11b80744cba015915e9db12ba324c01783817ecba153f62d20.
The new cubic-trace DQI phase kernel reduces each scan point from92to88Toffolis.
Full native102400shots PASS, zero value errors,380phase errors, Q961,
meanT1019865.8565, cost1598204229384.9265. Both archive audits and an independent
full source rebuild pass. Previous9.3101204888archive is preserved in
revisions/compact_product_slim/. Inverse and chemistry are unchanged.
SeeCUBIC_TRACE_DQI_RESEARCH.md. An87Tkernel is exact but has a worse native
finite-sample repetition cost and is excluded. All research jobs are terminal.

GROUPED_FACTOR_RESEARCH.md records a17.02mHa floating-point R28fit whose
quantized dictionaries remain invalid. PHASE_AND_GROUPED_RESEARCH.md records
new cubic arithmetic, Boolean and Walsh experiments. LOOKAHEAD_RESEARCH.md
records a pivot-robust quantum pair map and nonwinning direct-series blocks.
The following older observations are historical.
A further product
reduction of approximately3.38 across task costs is needed; this can come from
combined improvements. Existing local minima and rough resource estimates do
not prove the target impossible. The preceding9.2438728734 package is preserved
under `revisions/affine/`.

Latest selected change: **ANGLE_GAUGE_RESEARCH.md**. An exact signed Givens
gauge compresses lookup payloads from8 to7independent bits per word while
keeping the public eight-bit XOR contract. Full native chemistryPASS50880,
cost712929326.0737408, and all11oracles PASS2048random incoming targets each.
Independent rebuild and ZIP audit pass. The preceding9.2439070508 archive is
preserved underrevisions/jacobi/. Inverse precision screens generated accurate
new tables but lost on actual circuit cost; none entered the archive.

## Exact column rotations

`jacobi_l1.cpp` searches every smooth interval and sign breakpoint of each
two-column rotation, minimizing the combined squared L1 column norm exactly
along that plane. It preserves `A A^T`, where A stacks the native factor weights
and identity coefficient. `check_jacobi_l1.py` compares60 random/sparse examples
against32769-angle grids. All checks pass, with Gram error below1.8e-15.

Twenty16000-pair passes reduce the checkpoint's two-body normalization from
294.4025805889514 to294.39650157216454. The720-row Gram error is7.2e-13.
All NPZ arrays rebuilt bit-for-bit using the same seed and commands. Full native
chemistry verification passes50880shots, cost722581182.7790365 and
lambda_eff48.89470055602427. The selected candidate is `candidates/jacobi_l1`.

For this fixed combined Gram matrix, triangle inequality gives the lower bound
`lambda_h2 >= sum(abs(A A^T))/4 =275.14137959044`, only6.54% below the old value.
This limits gains available from exact refactorizations of **this Gram matrix**;
it does not limit a new dictionary or a changed physical approximation.
See `results/rotation_gram_lower_bound.json`.

`results/offset_other_branch_probe.json` checks a scalar path from the selected
identity offsets toward zero identity coefficients, including extrapolations.
The tested path worsens the effective normalization; this is not a global proof.
`results/augmented_gauge_bounds.json` records similar fixed-Gram bounds for
earlier715-vector fits. Those fits cannot cheaply beat the checkpoint by column
rotation alone at their current offsets and gap.

## LLL/Babai integer-angle jumps

`lattice_babai.py` converts the real local curvature metric into an integer
Cholesky lattice, reduces it with LLL, and proposes nearest-plane and reduced-
basis moves. It is a proposal method, not an exact closest-vector solver. Every
proposal is rescored with the original quadratic; chemistry acceptance uses the
actual nonlinear rank-constrained projected loss. Moves can be many angle words
away, unlike the previous ternary QUBO proposals.

The independent barrier example in `results/lattice_babai_check.json` has every
single-coordinate unit move uphill, while the lattice proposal (-10,10) reaches
the exhaustive optimum. fpylll0.6.4 and cysignals1.13.1 are installed locally
under `vendor/lattice_python`; their installation metadata retains licenses.
Use the torch environment Python, after sourcing env.sh.

`lll_coordinate.py` implements the physical search. The initial12-row search
accepted an18-angle move of up to2words. Searching all700rows with weak damping
accepted64-angle moves up to45words. The best actual-reference two-body error
fell from20.255337mHa to20.241422mHa. This is still invalid against20mHa.
Pairing nearby atoms was weaker, reaching20.252084mHa. An arbitrary128-coordinate
global block generated no improving proposal in its tested configuration.

The newer trust-region run limits each coordinate jump to8words, permits128
actual trials per iteration, and tests all700rows. Its output and log are:

    candidates/lll700_trust8_q8_c350.npz
    results/lll700_trust8_q8_c350.log

That100-step run is terminal: best20.214431mHa. A subsequent long continuation
was intentionally stopped to change the search schedule, not because of a
timeout. `lll_coordinate.py` now accepts `--moves-per-step` and tests every
candidate against the current nonlinear objective after each accepted move.
The current run accepts up to32moves per iteration,128trials,700rows,trust8,
padding.001,500iterations. Its atomic NPZ writer keeps the best actual physical
error, saves accepted intermediate bests and also saves the final iteration.

Current live run, confirmed process **26361**, unified-exec session **10923**:

    candidates/lll700_multi_q8_c350.npz
    results/lll700_multi_q8_c350.log

At step44 the best saved error was **20.144238mHa**, still above20mHa. Poll the
actual process/log before deciding whether it is terminal. Do not restart solely
because a wait times out. All earlier wide, paired, trust8-100 and trust8-long
processes are terminal. NPZ saving is now atomic, so it can be read safely.

Newer main-run observation: step145 best20.121882mHa, still live. A separate
paired escape search runs from the frozen seed inparameters/lll700_pair_escape_seed.npz
(SHA2cbc6bd9...1384384), session29690. It uses96single rows and24paired blocks,
padding.0001,maxstep32,128trials,8accepted moves/iteration,100steps. Its output
and log use the stem`lll700_pair_escape_q8_c350`. Both jobs are independent
classical research; neither has established a valid lower-cost chemistry fit.

Latest status supersedes the preceding observations: paired escape finished
after26steps at20.114304597mHa. The original main run was deliberately stopped
atstep172(best20.119658mHa), and the full700-row search moved to the paired best.
Current live process **28022**, session **31184**, output/log stem
`lll700_pair_polish_q8_c350`, usesmaxstep32,padding.0001,256trials,32moves per
iteration,300steps. It remains above20mHa at initialization and is not a valid
new submission. This switch used terminal evidence of a better branch, not an
observation timeout.

The successful large-word moves show the search is reaching beyond the previous
unit/ternary neighborhood. They do not establish a valid new chemistry fit.
`results/lll700_ridge_sweep.json` shows that lowering ridge all the way to zero
only reaches20.241306mHa for the prior wide-run dictionary, so that dictionary's
miss is not simply excessive ridge bias. These new error estimates use the full
reference spectrum, including negative eigenvalues; older positive-part fit logs
can differ by a few microhartree and must not be substituted for native results.

## Larger architectural routes still open

**BLOCKED_EUCLID_RESEARCH.md is the current detailed decoder handoff.** New
quantum generic steps, quotient blocks, full armed sequences and completion
padding pass value/inverse/composition/scratch/phase checks. A large classical
quadratic-padding probe passes34850cases/626265roots with zero artificial
collisions, skipping150degree-gap paths. The full padded quantum sequence is
tested atN22. Exact N70 forward profiles expose a substantial cost problem:
completion padding adds276120T/direction in the current three-pass schedule.
Do not mistake the generic278370T arithmetic schedule for the padded whole
decoder. No new DQI circuit has entered the archive. A fixed-degree Sylvester
adapter for H is derived but unbuilt, and generic zero-pivot assumptions require
an audit before integration.

The normalized reciprocal's optimized variant now passes512quantum shots at
length35,Q423/T11166 and length70,Q843/T32622. Its Möbius conjugation map passes
512shots including zero constant T. These are components, not task submissions.

The serial compressed Euclid estimate is a bound on a particular implementation,
not on all decoders. Its875-qubit first step is quantum-tested. A blocked Padé/
half-GCD update, or another in-place moment-to-locator/evaluator map, could remove
the present decoder's extra polynomial bank. See COMPACT_EUCLID_RESEARCH.md.

Another explored formulation contracts the classical runtime mask polynomial
modulo the quantum locator and pairs its remainder with retained syndromes,
removing explicit Chien/Forney scanning. Straight companion Horner is too
expensive; no fast low-workspace circuit has been built. Do not count this as an
implemented improvement or an impossibility result.

The existing inverse's27-bit arithmetic is partly a consequence of its repeated
partial-product rounding. A different arithmetic/phase representation may trade
that error budget and workspace differently. Preliminary counts for existing
cubic schedules are not a lower bound on new evaluation architectures.

An optional async question asks the user for a public code link or per-task
breakdown of16.6. The work continues independently of a reply.
