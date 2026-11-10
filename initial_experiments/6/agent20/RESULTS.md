> Current handoff: [REPORT.md](REPORT.md). The packaged best is **12.046515352810363**, milestones **10–12** are saved, and a newer fully checked decoder awaits packaging for a projected combined **12.11607972269253**. The older accepted-state statements below are historical.

# Score-20 research

ACCEPTED: `best/submission.zip`, SHA256
`df518f8d20a0b4d4ea72700280f689bc17d0685471e0d1ddb2cc998b1c60a52e`.
Verified local default score **10.128965274112238**; both changed tasks under
their predeclared independent salts give combined score **10.144446072819372**.
The first-crossed archive `milestones/10/submission.zip` is preserved unchanged.
Evidence: `results/accepted_H_improvement.json`, original metric audit
`results/audit_best_H_metric.log`. The goal of20 is NOT complete;
milestones11–20 remain outstanding.

Goal remains active: produce a verified submission scoring at least 20 and save
submission archives at every integer threshold 10 through 20.

Authoritative initial baseline: round-5 integration/wave3/submission.zip,
SHA256 b521c86dfdf830ada6e0f9cbfc88a5e58c762be3a2feabaa50c0479324652fa1.
Full-default local score 9.811059643318059. Costs: inverse 182956.5271606425,
DQI 1414354472410.0034, chemistry 696942846.9115301.

Initial direction: build and measure the previously numerical-only narrow
quadratic-square inverse architecture. Keep the full-precision linear coordinate;
round only the curvature coordinate. First establish an independently exhaustive
integer model for arithmetic with narrow output and actual Booth-row rounding.

All new writes stay inside agent20. Host verifiers and prior artifacts are
read-only. Full verifier runs use /tmp/r6-fullverify.lock. Circuit components
from rounds 2–5 retain their original attribution; see their source docstrings.

## Verified inverse improvement

The exact radix-four partial-product rewrite in `booth_adder.py` uses the
three-AND Karatsuba identity for adjacent coefficient pairs. It reuses the
cleared Booth difference wire and repays every internal measurement phase.
The numerical approximation and synthesis margin are unchanged.

`candidates/rows_best.kmb`, SHA256
e8a5ab6b8af80154bec09b457d01afc5ae01b3c1f35a15b021a51af828c5a831:
90 qubits, 1693 static Toffolis, 34 rotations. Both full-default native runs
pass 200540928 shots: default cost 173956.46294113316, independent-salt cost
173956.67613342477. Combined with unchanged round-5 pieces this scores
approximately9.9774. Evidence: `results/rows_*_verify.json`.

The independent exact-phase native diagnostic also checks 42334 distinct
inputs under two measurement seeds (84992 shots). An 89-qubit schedule has
1709 static Toffolis and full-default cost 173447.66648015822; its independent
salt also passes200540928 shots at cost173447.6937108191. The combined verified
score is9.98716898063121. This is selected by resources, not by measurement seeds.

## DQI retained inversion

Porting the prior retained-inversion primitive to the fixed accepted gap
schedule saves 5256 static Toffolis at unchanged Q1004 under a measured
retention limit of 980. The looser limit1004 reaches Q1016, and is rejected.
`candidates/retained_inverse_dqi.kmb` is built at 1167282 static Toffolis.
Full-default native verification REJECTS this port: value errors0/102400,
phase errors50501/102400. No cost or improvement is accepted. Opening and
selected-cofactor paths swap source and reciprocal; retained norm intermediates
then describe the old source. `safe_retained_inverse.py` preserves the original
primitive at those mutable sites and reserves each immutable body's scratch.
The corrected N14/M127 complete oracle passes4096 native shots with zero value
and phase errors (`results/safe_small_verify.json`). A new resource probe runs
before selecting the corrected full-size build.

The corrected source-preserving build saves9792 static Toffolis at Q1004:
core791596, complete1162746. Both full102400-shot native checks PASS:
default cost1396130633850.919, zero phase errors; independent cost
1402542281832.912, one phase error; both zero value errors.
Artifact `candidates/safe_retained_dqi.kmb`, SHA256
ab69453dae0ef73a42d1065e962baa0f050428565fd32e21e14c934241915318.
The smaller oracle and full verification validate the source-preservation fix.
Evidence: `results/dqi_safe_*_verify.json`, `results/accepted_milestone10.json`.

## New field arithmetic component

`cube_field.py` implements a 28-AND GF4096 cube using seven GF16 evaluation
points, four-AND GF16 cubes, exact interpolation and linear basis maps.
It passes all4096 field values with arbitrary target words under two native
measurement seeds (8192 shots); inputs, scratch and phase all restore.
Artifact: `candidates/cube28.kmb`; evidence `results/cube28_component.json`.
`cube_arithmetic.py` applies this to Frobenius-related products and considers
cube evaluation in series inversion. Full decoder resources and native
integration are still under investigation; no score claim for this component.

The serial C-fraction generic-chart resource profile improves from1450238 to
1376094 complete static Toffolis (878-qubit arithmetic core, before the scan
layout). Extraction falls503655→466583 per direction. This still lacks zero-head
chart repair and is not a competitive accepted decoder. The fixed accepted
gap core exposes only small cube instances; its measured cube gain is modest.

## Chemistry and further search

Exact-word rank-C capacity at B719 is confirmed below20mHa down toC340.
However, merely refitting and aligning to original physical generator subsets
gives lambda_eff about57.24 atC340, which would lose despite the smaller rank.
The initial PCA gauge is much worse (lambda_eff about103). These are scoped
weight/gauge tests, not exclusions of better dictionaries or constrained fits.

Quantization-in-the-loop training starts from the accepted words. Two initial
Adam controls retain no changed-word improvement: first2000 steps are unstable;
the smaller-step20000-step run briefly improves exact error2 to19.9375mHa
without changing words, but later word changes lose feasibility. All best states
and histories are saved. The next control projects the actual SOS gap in float64
at each checkpoint and constrains raw Lambda, rather than allowing the effective
normalization objective to exploit invalid gap values. That10000-step run is
active; it is a numerical research experiment, not an accepted submission.

The reoptimized fixed-arithmetic one/two/three-slot schedule saves only340
projected Toffolis at the same peak. No new full build is warranted from that
small resource-only change. `results/safe_schedule_frontier.json` preserves
the whole frontier; no accuracy claim is made for alternative schedules.

The rounded-square inverse model passes all3133452 inputs with fractional
precision1–4 (`results/square_models.json`), but no whole circuit was built.
The existing Horner architecture produced the first concrete improvement.

Milestone10 is saved and audited. The original host circuit/NPZ readers and
score function pass on `best/submission.zip`; every ZIP member is hash-bound
to complete native evidence. `package_best.py` saves future crossed thresholds
through20. `verified_inverse/submission.zip` preserves the intermediate9.9872
submission as a separately audited fallback.

## Subsequent blocked Jacobi research (not accepted)

The new Jacobi fraction f=a/(1+bx+x²g) consumes two moments per scalar
normalization. Its two-pair block shares one long reciprocal with an implicit
zero coefficient, saving substantial gates over blocked ordinary C fractions.
The compact head change c+=bd,d+=b consists of two reversible shears.
`compact_jacobi.py` emits Q915/core777004/whole1104604 static Toffolis.
Native component tests N6/14/34 pass exact forward/inverse value, phase and
scratch checks under two seeds. A complete N14/M127 oracle passes4096 shots
against its independent totalized model and restores all inputs/masks/scratch.
Forced unsupported charts give145 truth phase errors; supported charts have0.

An independent natural-law census finds121 wrong magnitude vectors/20000 for
both serial and blocked Jacobi. The expected random-mask phase errors and
host confidence penalty erase this generic architecture's resource saving.
This is a specific missing-zero-chart problem, not an accepted improvement.

`general_jacobi.py` prepares larger blocks' affine shift using a short series
B/A and a high-band middle product. This avoids the large expanded polynomial
cost. Block3 native N8/14/34 checks pass including phase and scratch. Full-size
resource profiles: block3 Q987/whole1037982; block4 Q1035/whole1027914.
Neither repairs zero pivots. No full-size generic Jacobi binary is accepted.

`variable_shift_inverse.py` is a new exact component for chart repair: the
quantum selector chooses r/(1+x²r) or r/(1+x³r), controlling writes only to the
extra final coefficient. Native N6/14/34 tests pass both selector values,
independent rational-series outputs, phase, clean scratch, and involution.
Its N34 cost is14237T at425Q including408 public series bits plus selector.

Chemistry projected training and robustness controls have finished without a
changed-word improvement. Fixed-weight integer descent with normalization
headroom lowers B703 two-body error to20.414969mHa, still above20mHa. This is
numerical progress only. Milestone10 and all accepted hashes remain unchanged.

The leading-zero chart is now emitted and checked. It reduces natural-model
wrong magnitude vectors to75/20000 versus161 serial on the same independent
stream. The missing second pivot remains. `shared_jacobi_synthesis.py` shares
ordinary/gap matrix work and repays early measurements after reconstructing
the original factorized head from low output polynomial coefficients. This
preserves common-factor information on armed charts. Native synthesis checks
pass forward/inverse, phase and scratch at N4/6/14/34 with deliberate gaps.

The current complete SMALL N14/M127 shared oracle passes4096 native shots with
zero value, scratch, model-phase, mask, or twice errors and zero truth phase
errors on that stream. SHA256
a567606ec6fa1d1ecd8a6763e2dc19ed6a52c14ae980b2b74ef11b95d397f59d.
Full-size resource profile is Q942/core984034/whole1311634. It remains
uncompetitive after the residual failure penalty, so no full build is selected.
The variable-shift component's V2 output aggregation costs10447T at436Q for
N34, replacing the earlier14237T/425Q version above.

The finite five-cycle chemistry alternation has finished at B703 error
20.144102913mHa, lambda_eff51, still above20mHa. The optimal-column exact-word
partition bootstrap also fails: R2/B351 and R2/B479 are around465mHa and237mHa
after PSD refitting. These are numerical research, not accepted chemistry.

The rounded-square inverse baseline is now fully built at Q122/static2063,
35 rotations. Native42334 distinct inputs/two seeds/84992 shots match the
independent exhaustive-passing integer model with clean input and scratch.
It loses to the selected Horner circuit, so no full native verifier is run.
SHA2565748a098be62abfd23f06f2892f6e1fce39b0c066245b92836e6e2d459928218.

New even/odd series inversion passes native exact reciprocal/phase/scratch/
involution checks at N2/3/6/8/14/34. N70 improves30158→27762T, but porting it
to the accepted DQI schedule saves only616T at the sameQ1004. No rebuild.

New in-place polynomial modular inversion uses continued-fraction reversal.
The improved unnormalized version computes a two-coefficient quotient's
entire tail action with one bilinear convolution. Native degrees1/2/3/7/17/35
pass1024 shots each, including forced zero pivots: exact independent totalized
model, phase, scratch and involution. The inverse-modulus identity holds on
the nonzero-pivot charts. After trimming an unused source coefficient,
degree35 costs73542T atQ892 including852 publicbits,
versus114138T/Q868 for the first monic prototype. This is a validated component,
not a full DQI improvement. Interior degree gaps and variable-degree/common
padding still need correct integration. Reports
`results/poly_inverse_bilinear_trimmed_checks.json`. Variable-prefix reversal
is now native checked through degree35: all padding degrees and arbitrary
paired inverse inputs/flags PASS; degree35 Q933/T78882 forward.

The same bilinear source-division rewrite improves in-place H preparation
from50040 to41220T. Native degree7/17/35 forward totalized-model and arbitrary
paired inverse checks pass. A new padding-aware first-gap repair has an
independent36000-fixture classical screen: all12 failures have repeated
locator roots or a zero locator root, outside physical decoder support.
The native repair emitter now passes degree7/17/35 tests including forced
gaps/all padding degrees and arbitrary paired-inverse inputs/flags. A binary
first-gap tag reduces degree35 width to957 at67956T forward. The independent
50000-vector decoder model has zero errors and170 repairs. Complete N14/M127
oracle native4096shots also passes with zero value/phase errors.

Full candidate `candidates/repaired_H_dqi.kmb`, SHA256
ac2d7a963c1b3e2f58f39a08000838296f4a6fd18d468cd35567351b211d41ed,
is built atQ1004/static1123942. Both full102400-shot native runs are now in
progress under the fixed default and predeclared independent salts. No new
score or replacement is accepted yet. See `CONTINUATION.md` for live status.

The verified score remains10.0304357867; milestones11–20 are still missing.
The score20 goal remains active.

## Accepted repaired-H full result

The completed repaired-H decoder passes both102400-shot native workloads:
default zero value/one phase error, cost1355783099065.6162; independent zero
value/phase errors, cost1349585429383.3171. It uses1004 qubits and1123942 static
Toffolis. The package passed the original archive reader and combined metric,
then replaced `best`; the old package is `verified_before_H`. This supersedes
any earlier in-progress wording in this chronological log. Reports and exact
archive digest are recorded at the top and in `CONTINUATION.md`.

## Finite multiplier decomposition screen

An exact GF(2) tensor flip/plus search completed59,395,257 transformations
(seed2026102501,120 seconds), preserving the full12x12x12 multiplication
tensor. Best rank remains42; no circuit change is selected. The search permits
six extra rank-one terms and restarts every10000 steps. This is a negative
result for that finite search, with no optimality claim. Source `tensor_flip.cpp`;
report `results/tensor_flip_search.json`.

## New normalization-aware four-block fit

A fresh R4/B240/C88 float32 run jointly optimizes the full signed tensor,
native Lambda and SOS gap, rather than rounding a prior high-cancellation
multiblock basin. Seed2026102503,12000 steps, learning rate.002 to1e-5,
checkpoints1000; `chem_frontier_fit.py`, tag `fresh_R4B240C88`. This is a
bounded numerical experiment, not a native acceptance or projected score.

The fresh R4/B240/C88 run completed12000 steps in170seconds: float error
63.5026mHa, rounded173.9854mHa, exact-word PSD-refit132.3009mHa,
refit lambda_eff78.6989. It fails; no circuit or score improvement is selected.
Report `results/fresh_R4B240C88.json`.

## Projective modular inverse

New source `projective_polynomial_inverse.py` removes alternating quotient
normalization by allowing Fnew=F/c^2, preserving locator roots while obtaining
Gnew=G^-1 mod F. Counter encoding and partial shifts exploit promised common
monomial padding; the explicit inverse covers arbitrary flags and operands.
Native all-degree and arbitrary paired-inverse checks are running for degrees
3,7,17,35 in `results/check_projective_poly_inverse.log`. No full-oracle or
score claim is made; interior degree-gap repair remains necessary.

The projective inverse checks completed: degree35 Q939/T71782 forward,
143564 paired, all padding degrees and arbitrary paired inverse PASS512 shots
per circuit. The accompanying in-place reciprocal passes exhaustive native
checks at66T (controlled91T); its norm chain is explicitly transported after
the source swap. Exact evidence is in `CONTINUATION.md` and the reports.

Sparse-pair fixed support screens fail with best residual lower bound247.30mHa.
A separate density/Schur structured bootstrap is running in
`chem_density_bootstrap.py`; no score claim is made.

The density/Schur bootstrap completed: its target truncation error is only
7.452mHa, but the selected exact-word dictionaries give591.083mHa after
10 PSD cycles, lambda_eff140.353. It fails. The subsequent conditioning-aware
B479→B511 growth also fails: best22.1378mHa, lambda_eff96.449; its lowest
normalization71.453 occurs at26.7321mHa. These are finite initializer/growth
screens, not shape-wide impossibility claims. No native chemistry submission
was built. All jobs completed; accepted best remains10.128965274112238.

The latest wave did not improve the accepted10.128965274112238 archive.
Native component tests now verify quotient-based Jacobi reconstruction with
all completion degrees through35 atQ933/T43255 per direction. Interior gaps
remain unsupported, so it is not a complete replacement decoder. Corrected
nested reciprocal-copy profiles lose on the accepted core; earlier outer-only
profiles did not exercise the copied inverse and must not be read as such.

A fresh chemistry representation now has352 positive semidefinite generator
matrices and7.36114mHa tensor error under the fixed native shift, constructed
by spectral/polar ADMM. Its intrinsic trace-based offset optimum is
lambda_eff47.91043 atgap4.08171. The compact R4/B240/C88 dictionaries tested
so far still fail: direct positive tensor fitting reaches38.083mHa continuous,
float64 L-BFGS-B35.950mHa, while rounded/refitted error remains much larger.
The tested best integer geometry has a92.355mHa union-span lower bound.
Finite native integer refinement of this new geometry is running; no native
oracle acceptance or score claim follows from the positive factor alone.
See CONTINUATION.md for reports, hashes, completed negative screens and
current process state. Goal20 remains active;11–20 are still missing.

## Latest mixed-quotient research

The accepted archive and score remain unchanged. Compact head frames,
high-half raw-head copies, and a seven-bit completion census now have native
component validation. Shared degree-one/degree-two quotient reconstruction
also passes independent rational identities, every completion degree, arbitrary
paired singular inputs, phase and scratch checks through nominal degree35.

The complete mixed SMALL oracle passes2048 shots, including293 arbitrary
syndromes, with exact model phase, mask/input restoration, scratch cleanup
and repeated application. The routed version uses Q280/T46432. Its168 true
decoding errors lie on unsupported charts; this is not benchmark acceptance.
Full-size emitted resources project to Q943/T967536, but a separate100000-
vector census finds658 wrong decoded vectors, so it is not promoted.

The original chemistry integer run finished at56.73mHa, above20mHa. Two new
factor ADMM controls and a finite global eigenvector atom-replacement variant
also failed to improve that result. A new three-slot cofactor parameter frame
passes native phase, scratch and inverse checks; its complete extractor is
unbuilt. CONTINUATION.md records exact scopes and results. All artifacts
remain local, and the goal of20 remains active.

The triple head and complete block are now built and native-checked. The
canonical head costs Q145/T2498 per direction; complete N70 extraction costs
425579T per direction and is too expensive. The generalized degree1/2
quotient reconstruction also passes native phase/inverse checks, including
arbitrary adjacent flags. A reversed-quotient cofactor frame and full b6
lookahead prefix block pass small native checks, but the N70 block profile
Q1000/T102376 loses to the existing prefix. None is benchmark acceptance.

Further weighted/global chemistry ADMM fails to improve56.729335mHa. A
finite R11/B64/C32 fit also fails; its saved exact-word dictionary has a
158.634857mHa unrestricted span error floor, excluding weight-only rescue
of those particular vectors. These are scoped negative results. Accepted
score remains10.128965274112238; only milestone10 exists. See the new
generalized-quotient research manifest and CONTINUATION.md for bindings.
