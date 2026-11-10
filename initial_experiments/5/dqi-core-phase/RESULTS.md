# Phase-aware decoding research

**Final outcome: no replacement promoted and no verified task-cost improvement.**
The accepted gap candidate is preserved read-only.
Baseline: Q=1004, mean T=1195802.1546875, R=1.1912604648115832,
cost=1435930717449.9397. No full-size 102,400-shot native acceptance run was
attempted for the uncompetitive approximate architectures below.

The strongest deliverables are the in-place C-fraction decoder, the in-place
Wronskian preparation and completed 46T inverse-only scan, and the two
phase-aware algebraic recurrences. Full proofs and implementation contracts
are in [DERIVATIONS.md](DERIVATIONS.md). Machine-readable whole-oracle screens
are in [results/resource_model.json](results/resource_model.json).

| Architecture | Full-size Q | Complete T | Status |
|---|---:|---:|---|
| Accepted gap incumbent | 1004 | 1,195,802.1547 mean | Existing full native verification, unchanged |
| C-fraction state, direct evaluator output | 884 estimated | 1,450,238 static estimated | Native complete small oracles; unsupported zero-head charts |
| Fixed-degree Wronskian state and 46T scan | 1004 assumed | 1,153,356 static estimated | Native full-degree components; composition/gap repair incomplete |

At the zero-error 102,400-shot confidence floor, the two estimates would give
only **1.06361x** and **1.03680x** improvements, respectively. These are
conditional budgets, not scores. Current independent classical model errors
erase both margins. The error-free screens deliberately do not pretend that
unimplemented gap repair has a known cost.

Final assessment:

- **A — Verified improvement:** none. No full native benchmark candidate is
  accepted, and no default/salted task score is claimed.
- **B — Validated research components:** a reversible C-fraction decoder uses
  the syndrome bank alone for coefficient extraction and B/Omega synthesis;
  native N6 and N14 full small oracles verify actual syndrome/runtime-mask
  contracts. A degree-35 H preparation costs 100,080T build/erase at 897Q,
  including 840 public bits; its nine-point native phase sandwich is
  918Q/100,494T. The 46T/66Q inverse-only body passes 61,440 shots. The Riccati
  and general-quotient Wronskian invariants are independently verified.
- **C — Scoped negative results:** serial reciprocal extraction is expensive,
  the implemented generic charts fail too often, a separately retained H bank
  loses the Q-squared objective, and the tested serial moment functional costs
  at least 11,833,500T before other work. These are implementation/storage
  results, not global lower bounds.
- **D — Unresolved hypotheses:** blocked C-fraction extraction with inexpensive
  degree-gap handling; a (B,H) locator update that avoids reconstructing A;
  and use of the Riccati tangent without carrying a second complete decoding
  state. None has an established low-cost complete circuit.

The completed gap, phase, and hybrid reports and relevant implementations were
read before developing recurrences. These establish storage/cost failures for
particular representations, not a global lower bound.

## Investigations

1. **Riccati tangent / evaluator coordinates.** For the generating series
   S(x)=sum S_i x^i, use the derivation D S=S^2 (mod x^N), i.e.
   D S_(2i)=S_i^2, D S_(2i+1)=0. On an exact monic Padé pair
   S=Omega/B (mod x^N), deg Omega<deg B<=N/2, the induced tangent is
   D B=Omega, D Omega=0. This identifies a specific derivative of the
   decoding map with the evaluator, without a final modular inverse.
   Reversible storage and arithmetic remain to be measured.

2. **Reciprocal-series / C-fraction coordinates.** Peel one formal coefficient
   a=f(0) and replace f by (a/f+1)/x. This stores one word per peeled moment,
   rather than two complete BM polynomial banks. The inverse local chart is
   f=a/(1+x*g). Generic zero leading coefficients require degree-gap charts.
   Rational synthesis and its reversible scratch must be charged.

3. **Moment adjoint phase.** Given locator B in the position variable, extend
   the moment sequence by its linear recurrence, pairing with a per-shot
   classical interpolation of masks. This eliminates evaluator inversion,
   but the charged recurrence can be much more expensive than a point scan.

Artifacts and reproduction commands will be recorded as they are verified.

## Verified checkpoint

The independent `check_algebra.py` run passes 384 RS cases in GF(16) and
GF(4096), with degrees 0..35, 74 degree-gap cases, 3,360 root-value checks,
10,080 runtime-mask checks, and 350 exact tangent checks. It also checks 300
Wronskian updates with arbitrary quotient degrees and 285 coprime Sylvester
factorizations with arbitrary degree gaps. `results/algebra.json` contains the
fixed independent salt and full degree-pattern census.

The useful Wronskian identity is **not claimed as a new discovery**: inspection
of the read-only round-4 `WRONSKIAN_RESEARCH.md`, `SERIES_RESEARCH.md`, and
`bezout_algebra.py` found the same target and a different 62,098T/direction
Sylvester implementation. The new recursive elimination here measures
50,040T/direction for the complete degree-35 H preparation, including the high
coefficient adapter: 100,080T for build and erase, 897Q including the 840 public
bits. This is an actual emitter/allocation profile, not a native full decoder.
It is not an apples-to-apples whole-task comparison with the earlier
E[0]=1/BM-template layout. Generic Euclidean charts remain a limitation.

`h_phase.py --point` passes 61,440 native shots and measures 46T/66Q. The
46T adaptation had been proposed in the round-4 report; the component is now
implemented and independently checked. The degree-7 H preparation/phase/inverse
sandwich passes 512 native shots, runtime-mask preservation, source and scratch
restoration, and repeated phase cancellation. The first rejected test omitted
the classical mask square map; its failure log is retained and the corrected
circuit passes. No rejected test is counted as a pass.

`cfraction_decoder.py` is a separate decoder architecture, not a transplant of
BM or compact cofactor updates. Its pair peel is

```
a = f0; b = f1/f0; T=(f-f0)/x
g = ((f1/T)+1)/x + b
f = a/(1 + b*x/(1+x*g)).
```

A normalized reciprocal transforms the existing tail in place; the two head
words remain in that same array. Recursive pair synthesis then converts the
C-fraction coefficients to constant-normalized B/Omega, also in place.
Zero pivots are armed to guarantee a reversible permutation; this does not
repair unsupported decoding charts.

The first complete small N6/M31 native oracle passes 1,024 restoration/model
shots including forced unsupported gaps: 47 chart exceptions and 23 truth
phase errors, zero truth phase errors on supported charts. All syndrome words,
private scratch, dynamic classical masks, and twice-executed phase restore.
The build-and-erase core measures 3,462T. The initial scan reserved separate
accumulators. The subsequent version reuses B1/B2/Omega0 via an invertible
linear evaluation map and passes the independent validations recorded below.

Exact full-size C-fraction arithmetic profile: 503,655T extraction,
57,664T synthesis, the same reverse costs, total 1,122,638T; 878Q core,
including all 840 input bits. These are emitter counts, not a complete
full-size accepted circuit. A runtime-mask Forney scan adds 80T/position.
The representation removes the additional locator/template banks, but repeated
reciprocal series are expensive. Generic zero-head failures must be counted.

## Final native component evidence

All counts below are component or small-instance native evidence, not the
full 102,400-shot benchmark. Runtime masks remain arbitrary classical inputs.

- `results/sylvester_t17_core-phase-sylvester-b.json`: 512 shots each,
  forward/inverse/both, complete coefficient outputs checked independently;
  885Q, 49,326T per direction. The separate forced-zero remainder-chain test
  passes 1,024 shots in both composition orders at t=3. Unsupported charts
  preserve data and scratch; correctness of the Sylvester formula on them is
  not claimed.
- `results/h_degree35_core-phase-H-b.json` and
  `results/h_degree35_core-phase-H-c.json`: 512 and 1,024 independent shots on
  the **same frozen bytes**, SHA256
  `9aa227c624ed8390f8cd79bb305c096fc91e67cca67fbf5648f67b86772eb5f4`.
  The one-point phase sandwich is 918Q/100,126T.
- `results/h_degree35_m9_core-phase-H-multipoint-b.json`: 1,024 shots,
  nine field points, including nontrivial beta^70 runtime-mask maps;
  918Q/100,494T. Source/mask/scratch restoration and twice-executed phase pass.
  Degree 6 and degree 7 tests separately cover even/odd polynomial layouts.
- `results/padding_pair_n35.json` and `results/padding_H_n35.json`: 512 native
  shots each direction/composition, arbitrary quantum completion flags and
  coefficient banks. Explicit nonmonic leading coefficients are included.
  Costs are 7,344T per direction for the A/B or H/B pair, 3,672T for H alone.
- `results/cfraction_n6_m31_core-phase-cfraction-b.json`: 2,048 complete small
  oracle shots, Q116/T5942. Zero syndrome, mask, scratch, model-phase or
  twice-executed restoration errors. There are 94 deliberately included
  unsupported charts and 51 truth-phase errors, all on those charts.
- `results/cfraction_n6_m63_core-phase-cfraction-c.json`: 2,048 shots,
  Q116/T8502, 95 unsupported charts and 53 truth-phase errors. All supported
  truth phases and every restoration/model-phase check pass.
- `results/cfraction_n14_m127_core-phase-cfraction-d.json`: 1,024 shots,
  Q212/T31434, 71 unsupported charts and 30 truth-phase errors. All supported
  truth phases and every restoration/model-phase check pass.
- `results/moment_k35_core-phase-recurrence-b.json`: 512 native shots,
  878Q/3166T phase sandwich (1,583T arithmetic per direction). It checks the
  direct moment overwrite, runtime mask and complete source/scratch restoration.
  The k=3 moment and Wronskian-transport components each pass 1,024 shots.

`results/cfraction_n6_profile.json` exactly reproduces the small native core's
Q/T before the same emitter-counter method profiles N70. The final algebra
run adds 40 explicit dual-series checks localizing the tangent to the first two
C-fraction heads, with 112 attempts and zero incorrect supported checks.

## Independent complete classical models and confidence screens

The final `census.cpp` streams use **squared** host eigenvector weights,
nonzero values, distinct uniformly sampled positions, independently seeded
masks, and no circuit hash or verifier sample stream. The first amplitude-law
censuses were rejected after finding the missing square; their parameter file
and JSONs are preserved under `rejected_amplitude_*`. They are not evidence
for the host-law conclusions. The executable now rejects weight files whose
sum differs from one.

| Fixed model stream | Shots | CF wrong pairs | CF mask phase errors | Padded H wrong vectors | H mask phase errors |
|---|---:|---:|---:|---:|---:|
| 20261024077 | 20,000 | 253 | 120 | 81 | 41 |
| 20261024078 | 20,000 | 291 | 140 | 78 | 40 |

All supported C-fraction pairs agree with the independent exact error-derived
pair. The H model constructs the correct compact A from exact RS data solely
as a validation oracle, applies the actual **armed** Sylvester transform,
then the literal inverse padding and cyclic unpadding networks. It checks
phase against independent masks. That model assumes exact core cofactors;
rare incumbent core errors are not included. Intermediate unpadded/padded
chart censuses are also preserved and are not pooled for confidence.

The host CP function is evaluated only at **explicitly hypothetical**
102,400-shot counts extrapolated from those model rates. These are sensitivity
calculations, not confidence statements about a native run. They price CF at
approximately **2.178–2.289 trillion**, and H at **1.762–1.769 trillion**,
versus the incumbent's 1.436 trillion. Neither route warrants a huge full-size
build or final native acceptance run in its current form.

At the current generic gate/width budgets, CF tolerates at most 25 hypothetical
phase failures in 102,400 shots before losing to the incumbent; H tolerates
only 12. Cheap generic-chart approximations therefore cannot be promoted.
Fixing exceptional charts at zero cost was used only as a favorable screen.

The H assembly charges four explicit pair adapters, 100,080T H lifetime,
10,790T inherited unpadding/preparation, the 824,740T core and all 4,095 point
bodies. Full routing/liveness is not integrated. Q1004 is a feasible target
assumption, and at Q1020 the error-free projected gain shrinks to 1.00453x.
This is why a gate-only result is not an accepted improvement.

## Obstructions and remaining opportunity

The **proved algebraic obstruction** is narrow: nonsquarefree B can make the
specific map A -> sqrt((AB)') noninjective. It does not rule out another state
encoding, extra retained information, or an approximate branch.

The **storage obstruction** is also narrow: retaining three literal degree-35
banks requires at least 1,260 quantum bits, before private workspace. Even an
optimistic zero-error screen with every other qubit free costs 2.145 trillion.
It does not establish an information-theoretic width bound.

The **expensive implementation** is serial C-fraction extraction: repeated
full reciprocals cost 503,655T forward despite fitting within 878 core qubits.
A 1.5x result at Q884, retaining the present scan/synthesis, needs extraction
<=292,699T forward; 2x requires <=164,158T. Those limits also have to accommodate
exception repair. Merely rearranging final phase preparation cannot meet the
400k–600k long-term goal while retaining the 824,740T incumbent core.

The strongest next hypothesis is a blocked coefficient-extraction algorithm
that retains the C-fraction/Riccati information in the syndrome array and
synthesizes B/Omega without materializing extra polynomial banks. The explicit
head-fold and tangent identities here supply its output and phase contracts;
a cheap tail update with complete gap handling is still missing. The separate
H quotient recurrence is useful only if the locator update can avoid an
inverse Sylvester reconstruction of A at every stage.

## Reproduction

Run in this directory. The environment disables bytecode writes into the
read-only completed experiments and directs temporary files here.

```bash
source env.sh
"$EXP_PYTHON" check_algebra.py
"$EXP_PYTHON" h_phase.py --point
"$EXP_PYTHON" sylvester_circuit.py --t 17 --shots 512 --salt core-phase-sylvester-b
"$EXP_PYTHON" sylvester_circuit.py --t 3 --restoration --shots 1024 --salt core-phase-restoration-a
"$EXP_PYTHON" h_phase.py --degree 6 --shots 1024 --salt core-phase-H-even-a
"$EXP_PYTHON" h_phase.py --degree 35 --shots 1024 --salt core-phase-H-c
"$EXP_PYTHON" h_phase.py --degree 35 --m 9 --shots 1024 --salt core-phase-H-multipoint-b
"$EXP_PYTHON" padding_adapter.py --kind pair
"$EXP_PYTHON" padding_adapter.py --kind H
"$EXP_PYTHON" cfraction_decoder.py --n 6 --m 63 --shots 2048 --salt core-phase-cfraction-c
"$EXP_PYTHON" cfraction_decoder.py --n 14 --m 127 --shots 1024 --salt core-phase-cfraction-d
"$EXP_PYTHON" cfraction_decoder.py --n 70 --profile
"$EXP_PYTHON" recurrence_components.py --kind transport --k 3
"$EXP_PYTHON" recurrence_components.py --kind moment --k 35 --shots 512 --salt core-phase-recurrence-b
g++ -O3 -std=c++17 census.cpp -o census
./census 20000 20261024077 > results/census_complete_a.json
./census 20000 20261024078 > results/census_complete_b.json
"$EXP_PYTHON" resource_model.py
"$EXP_PYTHON" audit.py
```

`parameters/host_weights.txt` is the unchanged host function's normalized
squared weight vector; `parameters/host_quadratics.txt` contains the legitimate
classical padding constants. `candidates/` preserves native KMB circuits and
construction tapes. `results/` preserves native checks, classical streams,
resource models and rejected diagnostic logs. No complete replacement artifact
is submitted. Final artifact hashes and dependency/baseline hashes are recorded
by `audit.py`.

## Attribution and integrity

The field emitter, multiplication/inverse recovery, optimized unit-series
reciprocal, isomorphism, controlled-transvection synthesis and Forney kernels
are imported read-only from completed rounds 2–4. The underlying disclosed
arithmetic is Apache-2.0 `andigu/uncompute`, pinned revision
`64b5a13aa249f1360846657392e6c66817b26d0f`; its existing source/license remain
under `../../3/vendor/uncompute`. The Wronskian target, a different generic
Sylvester factorization, and the proposed 46T body are acknowledged in the
completed round-4 Wronskian/series reports. The present recursive Sylvester
elimination, explicit-leading adapters, C-fraction pair extraction/synthesis,
reference checks, native wrappers, models and derivations are preserved here.

The verifier, scoring/confidence functions, reference inputs and accepted
candidate are read-only. No hash/sample-based circuit selection, salt selection,
confidence pooling, tolerance change or verifier modification is used. No other
task's submission or shared Git state is changed.
