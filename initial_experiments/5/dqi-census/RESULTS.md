# DQI decoder architecture census

**No replacement accepted. The >=2.5x objective is not achieved.** Eight distinct
architectures are priced in [PHASE1.md](PHASE1.md), including direct F', lazy
BM handoff, resident sigma/omega, non-EEA/BM Hankel elimination, wide Psi,
root lists, transposed residues and Wronskians. Full state/scan/reversal,
peak register budgets, T/R/objective assumptions and largest uncertainties
are given there. Numeric liveness/threshold errata are corrected explicitly.

Phase 2 produced **10 native component circuits,27,136 shots in their final
suites,zero value/scratch/phase failures**. Two fixed classical streams also
verify the F' identity at every degree0..35. The useful new component is a
low-memory quadratic reciprocal-product map,960T per overwrite with69 extra
qubits beyond its two input polynomials. Its degree35 extension is unresolved.

The phase-1 top targets were resident rational coordinates C and direct F' A.
Their straightforward native constructions lose; further low-width norm
constructions were tested rather than stopping at that first result. The
families remain open. No speculative resource target is presented as an
achieved decoder, and no classical failure count is substituted for host R.

The most plausible conditional target remains **C:Q940,T571921.3867**, cost
**602003162954 (2.6548x)** at zero native failures. This is an unimplemented
240k core+10k adapter budget, not an estimate assembled from working resident
updates. Its full coefficient/history representation is still the central
uncertainty. [BUILD_PLAN.md](BUILD_PLAN.md) gives the precise missing primitive,
slot ownership, resource gates and complete native acceptance sequence.
No design now has a component-supported whole-circuit >=2x projection, so
Phase3's implementation gate is unmet; the plan is explicitly conditional.

## Preserved selections and constraints

Reference DQI:Q961,mean T1019865.856484375,R1.6968460732519504,
cost1598204229384.9265;0 value/380 phase failures in102400 native shots.
Selected artifact remains read-only:
`../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb`, SHA256
`4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0`.
Combined reference9.3677; inverse182956.5272; current chemistry708156073.7876
banked by chem-audit, per user. Earlier reports'708547388.9553 chemistry value
is historical. No archive, chemistry, inverse, verifier, confidence or scoring
file was modified. Every new file is in dqi-census.

At most two single-thread component processes ran concurrently; BLAS/OMP/MKL
are capped at one. No full verifier ran. A future full run must hold
`flock /tmp/r6-fullverify.lock` and use at most two workers. There is no full
native DQI candidate in this directory and no default/salt acceptance report.

Read first: round-5 repair/phase/hybrid/core/dispatch RESULTS; round-4 REPORT
and Codex continuation; [Khattar et al.](https://arxiv.org/html/2510.10967v1)
and the [discovered decoder README](https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md).
The paper's in-register state sharing motivates the census; its leading-order
costs do not certify this finite-size objective. The public circuit optimizes
TQ with a fixed mask; this host charges RQ²T and arbitrary independent masks.
Only dispatch's RESULTS was consulted; none of its changing implementation
files was imported or built. It currently reports a correct classical transfer
model and a losing matrix/Bennett implementation, not a selected replacement.

## Phase1 ranking, with assumptions exposed

R0=1.1912604648115832 is the exact unmodified host CP result for ZERO failures
in102400 shots. It is conditional on a new circuit's actual native outcome.
Rref=1.6968460732519504 is a separate sensitivity scenario. Neither is a
bound on an unbuilt architecture's new failures. Expected active points after
70-point interpolation:4024.017333984375. All branches are charged by emitted
arithmetic, not by how often their input pattern appears.

| Architecture | Q budget | T budget/estimate | Conditional R0 objective | Ratio vs incumbent | Main unresolved item |
|---|---:|---:|---:|---:|---|
| C resident sigma/omega, planning target | 961 | 611921 | 0.6732e12 | 2.374x | alias-safe update and template recovery |
| A direct F',100k converter/direction allowance | 961 | 903858 | 0.9944e12 | 1.607x | genuinely in-place degree35 map |
| H Wronskian, generic existing adapter | 961 | 973383 | 1.0709e12 | 1.492x | total adapter and H division |
| B early BM/compact handoff,k1,zero-cost handoff | >=1037 | >=1012335 | >=1.2968e12 | <=1.232x | reversible representation bridge |
| E wide BM/Psi | 1837 | 442397.4 | 1.7784e12 | 0.899x | remove width with compatible lifetimes |
| F root list/indexed masks | 2200 | >=443290 | >=2.5559e12 | <=0.625x | canonical factorization; no syndrome copy |
| D general Hankel displacement layout | 2200 | 655241 | 3.7779e12 | 0.423x | compressed rank-revealing generators |
| G streaming residue | 1360 | >=12689102 | >=27.9586e12 | <=0.057x | fast low-memory transposed reduction |

These mix measured parts with explicitly identified allowances; they rank
conditional design targets, not posterior expectations or complete-circuit
measurements. A credible scalar expected cost for an algorithm whose central
map is unimplemented is unavailable. Machine-readable estimates and Rref
sensitivities are in [results/budget.json](results/budget.json).

Two stretch targets are kept separate:

| Unimplemented target | Q | T | R0 cost | Ratio |
|---|---:|---:|---:|---:|
| C core240k+adapter10k+80T scan | 940 | 571921.3867 | 0.602003e12 | 2.6548x |
| A trimmed core624206,counter layout,converter30k/direction | 948 | 738262.1907 | 0.790377e12 | 2.0221x |

A's additional gap repair G is unpriced in this stretch target. That is a
condition for further work, not a proposed zero-cost correction. C's native
reference does not attain its target. At C's940Q/T target,15 failures give
R1.2425853,cost0.627940e12;24 give0.638647e12;25 exceed the2.5x allowance.
At380 failures the cost rises to0.857501e12 (1.864x).

## Phase2 C: native rational-coordinate update

The new reference exactly implements

```
V <- V / sigma_old mod x^n
sigma <- sigma + alpha*x*tau
V <- sigma_new * V mod x^n
```

sigma0=1. The whole operation is an involution at fixed alpha,tau. Forward
checks use independent host-basis polynomial arithmetic. Roundtrip checks
cover arbitrary input polynomials, zero alpha, zero tau and zero V. All
sources, updated targets, physical scratch and measurement phases are checked.

| Native component | Q | T | Shots | Result |
|---|---:|---:|---:|---|
| n8/d3 forward | 195 | 1470 | 512 | PASS |
| n70/d35 forward | 1707 | 31962 | 256 | PASS |
| n34/d17 forward+inverse | 843 | 25116 | 512 | PASS |

**Widths are component totals.** n34 omits36 untouched input words; embedding
it in the actual70-word task already needs>=1275Q before decoder history.
The degree35 peak partitions exactly:

| Co-live register | Q |
|---|---:|
| V,70 words | 840 |
| sigma including explicit constant word,36 words | 432 |
| tau,35 words | 420 |
| alpha | 12 |
| emitter pivots | 3 |
| Total | 1707 |

Making sigma0 implicit saves12Q; it does not remove the sigma/tau duplication.
T31962=2*15246 for triangular transforms+35*42 for the correction. Even
70 uniform reference updates and their inverses, before producing corrections,
would cost4474680T; adding the scan and10k adapters gives4806601T,
R0 objective16.68445e12 at1707Q. This prices the tested uniform construction,
not a lower bound on every rational decoder.

**Resolved against the reference realization:** it cannot deliver the240–280k
resident core or the940–961Q budget. **Unresolved for C:** an alias-safe direct
prefix transition and a short correction descriptor. At a generic2k boundary,
k sigma+k omega coefficients already fill the processed2k words; a separate
full tau has nowhere to live. A zero logical row that physically stores sigma
must never be read as zero by a disjoint-register triangular routine. The
native component supplies no erasure proof for such aliasing. Consequently
C remains an open architecture, not a dropped candidate or a completed design.

## Phase2 A: direct F' and the charged conversion

At compact roots beta, e=beta^70/(A(beta)*B'(beta)). Set
F'=(A B')^-1 mod B and absorb beta^70 into the classical trace-dual mask.
The scan then costs11T per active point. The new independent identity check
passes two predeclared seeds,288 cases,all degrees0..35,5040 actual roots and
1179360 independent mask words. It also checks f(f(A))=A on the promised
split-squarefree/unit domain. This is classical algebra evidence, not decoding
or circuit acceptance.

First concrete route: optimized Frobenius inversion in the quotient algebra.
All modular squares include quantum-modulus reduction; none is called free.

| Native component | Q | T | Shots | Result |
|---|---:|---:|---:|---|
| degree35 modular square XOR | 2091 | 17556 | 256 | PASS |
| degree3 F' overwrite | 999 | 55440 | 512 | PASS |
| degree3 arbitrary-state paired inverse | 999 | 110880 | 512 | PASS |

The exact degree35 deterministic emitter profile is11751Q/1799952T per
converter direction, or3599904T lifetime. It uses11 charged squares+10 charged
products per inverse evaluation, explicit compute/copy/uncompute, and two
reciprocal-product XOR evaluations per overwrite. A square costs17556T;
product23352T; inverse-XOR853272T. The decomposition exactly agrees with the
native small-degree circuits and degree35 primitive. This improves the prior
schoolbook arithmetic screen but remains uncompetitive.

| Co-live ownership at Frobenius peak | Q |
|---|---:|
| Existing A and B | 840 |
| Retained output/work bank | 420 |
| Derivative and A*B' | 840 |
| 21 inversion intermediates | 8820 |
| Convolution bank | 828 |
| Emitter pivots | 3 |
| Total | 11751 |

Even omitting decoder flags/degree routing/gap repair, the complete projection
is4303762T at11751Q, cost707.9534e12. **This implementation is rejected.**
The prior phase experiment already showed the broad Frobenius route's width
problem; the specialized square probe supplies new exact reduced counts and
does not pretend that this was a previously untested architecture.

A second construction avoids Frobenius entirely at degree2:

```
B=x²+b1*x+b0, A=a0+a1*x
u=a1*b1
N=a0²+a0*u+a1²*b0
F'=(a0+u+a1*x)/(b1*N)
```

It uses one scalar inversion per XOR evaluator. XOR/swap/XOR overwrites A;
the work bank remains allocated until reversal even on singular states.
Both versions test promised split B/unit A on forward shots and arbitrary
B/A/work on paired inverse. Controls include zero polynomial, nonunits,
repeated roots, and all4096 b1 values. Singular forward values are totalized
as a function, not asserted to be modular inverses.

| Native norm construction | Q | Additional Q beyond existing A/B | One-direction T | Forward / paired-inverse shots |
|---|---:|---:|---:|---:|
| Retain norm | 129 | 81 | 792 | 4096 / 4096 |
| Release/recompute norm, virtual square operand | 117 | 69 | 960 | 8192 / 8192 |

All pass, including scratch and phase. At the low-width peak:existing A/B48Q,
retained old-A work24Q,u+denominator24Q,scalar inverse pool18Q,emitter pivots3Q
=117Q. The entire69Q extra allocation is charged;45Q scratch alone is not
the conversion overhead. This fits the74Q target ONLY AT DEGREE2.

**Resolved against whole-bank XOR/swap/XOR at degree35:** its retained work
alone adds420Q, so even a free arithmetic evaluator cannot meet74 extraQ.
**Still open:** a genuinely in-place modular inverse with smaller blocked
workspace. The norm base case does not provide a quadratic tower for an
arbitrary odd-degree35 polynomial. Root factoring/order must be paid if used.
No17.5-times extrapolation of the degree2 count is accepted.

Break-even corrections: without interpolation,77*4095/2=157657.5T/direction.
If both scans interpolate,77*4024.0173/2=154924.6674. Against the actual old
mean88T scan, the interpolated11T scan permits158003.9154T/direction just
to reduce T. Those thresholds are NOT the2x budget: at unchanged961Q/core,
only11248T/direction permits2x even at R0; C=0 gives just2.0639x. A2.5x
requires further core/width changes. budget.json includes exact limits.

## Other candidates: sharpened checks and remaining uncertainties

**Lazy BM handoff B.** Inherited exact phase/liveness records supply a new
conditional splice screen at k1,4,7,10,13,16,19,22. The prefix is already lazily
allocated. No new BM/compact bridge was built. At k10, prefix38705T+compact
suffix477574T+scan370152T=886431T BEFORE handoff, peak>=1145Q. Its optimistic
R0 objective is1.384403e12 (1.154x). The cheapest screened zero-cost handoff
is k1:1.296846e12 (1.232x). Thus a later handoff's T saving does not beat its
Q² charge in these frames. This rejects those compatible-state,zero-cost
splice budgets as routes to2x; it does not resolve whether a different
representation can have a cheaper bridge and smaller peak.

The initial table's1140Q was instantaneous round20, not the cumulative prefix
peak. The correct r0..19 maximum1145Q partitions as840 input,144 locator,
36 private template,19 stack,64 fixed/anonymous,42 other temporary. The
full35-word locator storage floor1272Q applies only to the old representation.

**Wide BM/Psi E.** Read-only pinned liveness,first1837Q peak at T175139 in
`euclid/R`, partitions exactly:840 syndrome,444 template,432 locator,48 replay
pool,19 stack,12 EXT,12 one-word,12 discrepancy,8 counter,6 event record,
1 preswap,3 pivots. The840 syndrome words already contain omega/cofactor
workspace. They are not additional free tail qubits. Q1837,T442397 still
loses even at host R0. The scan's23T count is not a substitute for a complete
low-width lifetime.

**Hankel/PGZ D.** General displacement rank2 needs two generators on each
side; symmetry in characteristic two does not remove diagonal data for free.
The explicit2200Q realization loses even with optimistic quadratic arithmetic.
A compressed rank-revealing representation remains unresolved; plain dense
PGZ needs14700Q for35² field entries alone. Singular/rank-deficient matrices
are routine at weight<35 and cannot be put in a negligible failure bucket.

**Root-list F, residue G, Wronskian H.** Their explicit losing layouts are
priced in PHASE1 and budget.json. Quantum-addressed arbitrary classical masks
still need lookup gates; quantum divisors still need reduction; the generic
Wronskian adapter still needs total gap handling and H division. None of these
missing operations is waived by the rare occurrence of an input or by the
fact that its answer is mathematically determined.

The public18.9 top is evidence that the product of task ratios can improve.
It neither identifies which task improved nor proves that these DQI targets
have a realizable decoder. Failures here are evidence about specific
representations and circuits, not an architecture-wide optimum claim.

## Reproduction and artifacts

Run from repository root; env.sh disables bytecode writes into old sources.
The commands below are component tests, not shortened host acceptance tests.

```bash
source initial_experiments/5/dqi-census/env.sh
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py profiles
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py rational --n 8 --d 3 --shots 512
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py rational --n 70 --d 35 --shots 256
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py rational --n 34 --d 17 --shots 512 --repeat 2
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py square --n 0 --d 35 --shots 256
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py fprime --n 0 --d 3 --shots 512
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py fprime --n 0 --d 3 --shots 512 --repeat 2
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py quadratic --n 0 --d 2 --shots 4096
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py quadratic --n 0 --d 2 --shots 4096 --repeat 2
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py quadratic_low --n 0 --d 2 --shots 8192
"$EXP_PYTHON" initial_experiments/5/dqi-census/native_components.py quadratic_low --n 0 --d 2 --shots 8192 --repeat 2
"$EXP_PYTHON" initial_experiments/5/dqi-census/fprime_identity.py
"$EXP_PYTHON" initial_experiments/5/dqi-census/splice_screen.py
"$EXP_PYTHON" initial_experiments/5/dqi-census/budget.py
"$EXP_PYTHON" initial_experiments/5/dqi-census/audit.py
```

candidates/ contains only native COMPONENT KMBs. results/ contains their exact
counts and checks, budgets and the first numeric checkpoint. PHASE1 predates
component construction; later numeric errata are marked here. BUILD_PLAN is
conditional R&D guidance, not an integration readiness claim. [ATTRIBUTION.md](ATTRIBUTION.md)
identifies Apache-2.0 reuse and local/original work; source and artifact hashes
are retained. No full decoder was built, verified, promoted or selected.
