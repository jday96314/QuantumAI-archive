# Agent 14 research checkpoint

Consolidated account of successful and unsuccessful approaches: [Optimization report](OPTIMIZATION_REPORT.md).

## Goal achieved — audited local score 14.064304039685362

The complete `submission.zip` passes the archive audit and is frozen in
`milestones/14` and `checkpoints/14.064304039685362`. Its SHA256 is
dbcc81d3dcd80e416f203ff63f3539db65347cdfd27636259f68192224ff3396. The score uses the unchanged host formula
and full default-shot native reports bound to the exact extracted archive.

| Task | Selected artifact | Q | Full-shot cost | Verified shots |
|---|---|---:|---:|---:|
| Inverse square root | pair_curve_live_q85_t1763.kmb |85|170242.24959913077|200540928|
| DQI | raw_cf_destructive_n70_m4095.kmb |892|515978772816.4064|102400|
| Chemistry | round5 chem-census structured, unchanged |454|696942846.9115301|50880|

Inverse maxT1763, meanT1751.9997843731928,34 rotations; the original
approximation and phase margin are unchanged. Paired Booth convolution saves
ANDs; an exact52-product decoder reconstructs curvature from live C0 low21,
coarse address8, and slope high4, allowing the full curvature register to be
freed during the linear multiplication. The small preflight and all full
prescribed shots pass. Candidate SHA:
0bfbd90e829818d6e4da838c9bf08af50a199710e076efa0b5b6cd988269742e.

DQI maxT534575, meanT534568.8609277344. The new quotient consumes its
remainder evaluation and reuses those24 wires for norm/inverse work. Its
measurement phase is repaid by the bilinear identity num=head*den, without
recreating num. Restoring the numerator in the evaluation basis avoids that
extra multiplication. An explicit42T correction accounts for the two absent
proper-outer coefficients. A tower-inverse lifetime change reuses the norm
cube workspace. The quotient costs348T with16 private wires, versus664T and
20 private wires previously. Full synthesis native256 PASS, Q878/892,
T34750/45376; small complete8192 PASS with zero errors. The full candidate
saves12572T and4Q. Candidate SHA:
9e36454e027e14963dee25ed7b2595b2ab752df3d6982dd231ca32d93eaf8fd6.

Default DQI102400 PASS5value0phase. A separately declared fixed independent
sample also PASS5value0phase, cost515978842691.4708, salt
agent14-raw-cf-destructive-independent-a. These are separate samples; the
archive score uses default. The8 saved regression failures match the prior
circuit exactly under both fixed measurement seeds. Errors are included in
the unchanged verifier confidence/repetition penalty. No hash/salt tuning.

Evidence: `results/package_goal14.log`, `results/archive_audit.json`,
`results/pair_curve_live_q85_t1763_verify_default.json`,
`results/raw_cf_destructive_n70_verify_{default,independent}.json`,
`results/all_milestones_final_audit.json`. All milestone ZIP hashes and first
crossing scores10–14 were checked; earlier milestones were not overwritten.

| Milestone | First complete score |
|---|---:|
|10|10.268368324127126|
|11|11.02830129721277|
|12|12.078134489004718|
|13|13.202469964520382|
|14|14.064304039685362|


## Prior checkpoint: 13.827937894509976

`submission.zip` is audited and frozen under
`checkpoints/13.827937894509976`. Archive SHA:
5f6617bf09cba61db70096c9559bc8de9ac26845e6a15ff031717ba54de6cb94. First milestones10–13 remain unchanged;
only14 remains. The13.584101391938912 and earlier checkpoints are preserved.

Inverse now uses `candidates/pair_booth_carry_q89_p17.kmb`, SHA
ccdb5e088e7c94119e394106b3e52b45415f77cb422a5e8ad4b1accb66626efb.
Full200540928 PASS, Q89, maxT1709, meanT1697.999656399316,
34rotations, score173447.63819059133. The approximation is unchanged.
`booth_pair_carry_adder.py` replaces four adjacent convolution products by
three, includes rounded carry-boundary products, and borrows a zero Booth
parity wire for the product/adder workspace. Native kernel checks pass.
The earlier paired variant Q89,T1727 also passed all200540928shots at
175049.68185426682; its complete candidate is retained but was not packaged.
DQI remains threejetQ896,T547147,cost532860186278.65674; chemistry unchanged.
Evidence: results/package_pair_booth_carry.log and results/archive_audit.json.

Current research: exact quadratic curvature replay from live coefficient
registers. `curve_register_restore.py` kernel passes512 table cases, but
`pair_curve_invsqrt.py` failed full-domain preflight because delayed C0 high
bits were left live while borrowed as product workspace. No failing circuit
was promoted. `curve_live_restore.py` instead uses C0 low21, address coarse8,
and the live slope high4, with a53/52-product exact decoder fitted on all
possible live slope values. This preserves original delayed-bit timing.
New `pair_curve_live_q85_t1763.kmb` preflight is running.

Cubic Horner exploratory fit f6 widths27/19/10/5 passes exhaustive host-law
numeric fitting with167 leaves, but its complete prototype Q92,T1982 loses
and is not selected. Narrower tested cubic families failed numerical accuracy.
Radix8 initially failed because binary magnitude3 partial products require
integer carries; the direct XOR product combination was invalid. It is not
selected. Chemistry B671 tree/exchange jobs remain terminal without progress.


The goal is achieved. All first-crossing archives10–14 are frozen.

## Previous audited selection

`submission.zip` scores **13.485543697862889**, audited and frozen in
`checkpoints/13.485543697862889`. First milestone13 stays at13.202469964520382.

| Task | Selected component | Full-shot cost |
|---|---|---:|
| Inverse square root | round4 zero_head_affine, unchanged |182956.5271606425|
| DQI | candidates/raw_cf_lowtemp_n70_m4095.kmb |544628838530.43414|
| Chemistry | round5 chem-census structured, unchanged |696942846.9115301|

DQI Q896,T557645,SHA7cc5863d695ab940fd9f4190367b4c3b09fa8ac718ce7a848d56a35be708ada2.
Default102400 PASS6value0phase,cost544628838530.43414; fixed independent
PASS7value0phase,cost546103877154.05615,saltagent14-raw-cf-lowtemp-independent-a.
The8capturedseedlessdefaultfailures matcholdraw andlowtempexactlyunder2fixed
measurement seeds. Separate samples; packageusesdefault. No hash/salt tuning.
ArchiveSHA30a83c55a948e1d4f045d9bb6aa1bf812d17b85d8388a751570a59dc557c71b8.
Evidence:results/package_raw_cf_lowtemp.log,results/archive_audit.json,
results/raw_cf_lowtemp_n70_verify_{default,independent}.json.
Only milestone14 remains; all earlier thresholds and checkpoints are preserved.

## Prior DQI checkpoint: three-field jet

The selected complete archive now scores **13.584101391938912** and is
frozen under `checkpoints/13.584101391938912`. First crossings10–13 remain
unchanged; only14 remains. The previous13.485543697862889 archive is retained.

`raw_three_jet_model.py` replaces the four-field infinity jet by (p,q,d)
and chooses endpoint tail N=[d/q,p], B=[q], with the physical q=0 chart
N=[p], B=[] handled reversibly. The synthesized final denominator is already
monic, removing both70-field normalization passes. The monic outer fold
also avoids leading-field inversions. New synthesis N70 native256 PASS,
Q878/896, T34750/57948. Endpoint512 PASS and each N6/N14/N70 model2048 PASS.
The small complete decoder passes8192 shots with zero value/phase errors.

Full DQI Q896, T547147, meanT547140.3999707032, candidate
`candidates/raw_cf_threejet_n70_m4095.kmb`, SHA
47a2683ab44a921fafaaee5c17245d3454ec8311f5cddcefccaa460f1c2d1200.
Default102400 PASS5value0phase,cost532860186278.65674. Fixed independent
102400 PASS5value2phase,cost535823510019.331, salt
agent14-raw-cf-threejet-independent-a. The initial independent launch used an
unsupported environment variable and was stopped; only the corrected
`--salt` run is independent evidence. All8 saved common failure cases match
oldraw/lowtemp exactly under both measurement seeds142103/142104.
Inverse and chemistry remain unchanged. Exact archive SHA:
cdcee9af05716bd027065f80b26c9a80ec9f1f2bd000388364c8c662fe1a10ee.
Evidence: results/package_raw_cf_threejet.log, results/archive_audit.json,
results/raw_cf_threejet_n70_verify_{default,independent}.json.

Chemistry: B671 fixed-size15x8 dictionary exchanges and all88 native-word
NNI tree neighbors failed to improve21.10577661888594mHa. These jobs are
terminal and invalid fits were never promoted.

## Latest verified kernels and rejected alternatives

* `cubic_mask_phase.py`: classical Gaussian elimination factors multiplication
  by the runtime trace-dual mask, then applies its recorded row operations to
  scale the quantum divisor in place. A cubic GF16 trace-inverse kernel costs
  31T/Q31. All4096sources against12mask bases, zero and mixed masks PASS57344.
  `cubic_low_workspace_inverse.py` combines existing66T creation and31T mask
  repayment:97T,Q44 in the8192-shot native inverse-copy test.
* `cubic_mask_phase_joint.py` additionally reuses the round4 frozen16-AND joint
  quadratic cover (copied read-only into localparameters). Native57344 PASS,
  Q31,T30. It is not yet used in a complete candidate.
* `lean_emitter.py` builds independent parity operands inside source registers,
  reverses both basis changes, and retains only one clean product fanout bit.
  Native large synthesis/conversion checks pass. This saves two emitter wires.
* Reciprocal reverse extraction frame saves60T/row. One-product repair
  lookahead tests current actual b=0 and currentrawb!=nextrawb, eliminating a
  second multiply. `schur_extract_cubic` nativeN6/N70 PASS,Q899 atN70,
  T155004/177229 (read latest buildJSON for authoritative exact counts).
  Extension inversion uses192T forward/223T reverse,Q61;49152nativePASS.
* `schur_dynamic_division_sequential.py` uses the low quotient output as norm
  workspace, then recovers that output through a selected nonzero divisor.
  Native16384 including bothaxes,den0 and quantum nodeflags PASS. It trades
  gates for12scratchqubits. Outer synthesis now holds25head/completion ghost
  bits and pays them on reconstruction, removing the forward outer phase.
* Counter-free completion initially confused a real monomial gcd with a fake
  marker. Its full candidate passed but lost at765505462049.9541 (14value,
  5phase). The same19captured cases distinguish this bug from the reciprocal
  and quotient changes. The fix also checks the preceding polynomial's
  leading ratio; all14reconstructable common cases now match the old converter.
  All-padding nativeN3/N7/N35 tests pass. The corrected lean complete candidate
  passes default8value0phase and independent9value0phase, and was audited at
  12.148365926523397 before the current cubic selection.
* The exact chemistry split which replicates U719 while partitioning ensemble
  columns has invariant modeledT and greaterQ; no advantage. R4B255C90 noise
  training fails (.0374023continuous,.100659refit). New angle-lattice ADMM
  reduces that refit error to .0691379, still invalid. Unscaled ADMM was poorly
  conditioned; angle-scale40 and longer continuation were tested. No such
  invalid fit is promoted. A new R2B383C200 continuous fit is running.
* Quadratic-cover searches on GF4096 multiplication retain42ANDs, including
  cap2/4/6 walks. A joint GF16 front-end search retains16 so far. No new field
  multiplication saving is claimed.

## New raw continued-fraction architecture — full run accepted at milestone13

`raw_cf_model.py` removes scalar normalization and node swaps entirely. Raw
heads li=a_i+b_i*x obey R_i=li+q_i/R_(i+1). A projective two-term infinity jet
P0,P1,Q0,Q1, initially0,0,1,0, determines the terminal affine tail. Its update
is P0'=Q0,P1'=Q0+Q1,Q0'=P0+bQ0,Q1'=P1+bQ1+aQ0. A finite terminal can be
represented without division as N0=P1Q0+P0Q1,N1=P0Q0,B0=Q0²; a completed
infinite terminal uses N0=1,N1=B0=0. The numerator/denominator recurrence is
N<-(a+b*x)N+qB,B<-oldN. Its product is a degree-one convolution, costing one
field product per output coefficient rather than two per input coefficient.

Algebra checks PASS2048cases atN6,N14 and4096atN70, including1933forced scalar
pivots atN70. All sampled padding is squarefree, involving irreducible node
quadratics, and B has degree34 or35. A conditional x factor makes it degree35;
its zero root is outside the4095-point phase scan. Node-quotient recovery of
allheads passes512cases eachN6/N14/N70, including4408zero divisors atN70.

Native components: `schur_extract_raw` N70Q877,T114240/132685,512PASS;
`raw_cf_endpoint` N70Q926,T4547/8842,512PASS; `raw_cf_step` m36Q901/926,
T1554/2105,512PASS; `raw_cf_synthesize` N70Q926,T37848/61428,256PASS.
`raw_cf_normalize` pads even degree and normalizes coefficients while keeping
its scalar for the inverse. The generic squarefree converter is
`derivative_inverse_raw.py`; the middle uses the existing11T phase scan.

Complete `candidates/raw_cf_n6_m63.kmb` PASS8192,zeroerrors,Q158,maxT7732.
`build_raw_cf_oracle.py` is building and verifying the complete N70/M4095
candidate. No large-instance score or milestone is claimed until that finishes.

## Implemented arithmetic improvement

`triangular_dp.py` combines the earlier zero-row/column trimming with a new
search over every recursive split of each in-place unit Toeplitz map. The
algorithm selects the minimum number of actual middle-product field products;
ties favor balanced recursion. Prefix-preserving products and Frobenius powers
are supported. The decoder's previously selected gap-aware stage schedule is
unchanged.

`check_dp_arithmetic.py` checks four varied operand shapes, including nonzero
prefix boundaries, known-zero holes, and powers 1, 2, and 4. All 16 native
forward/inverse/composition circuits pass, with arbitrary operands, restored
scratch, and exact phase. `build_dp_decoder.py` applies the arithmetic to both
outer products and directly constructed local-head arithmetic.

## New field-cubing primitive

`gold_interpolation.py` constructs an exact **28-AND GF(4096) cube**, compared
with the generic 42-AND multiplication. Represent an element as a quadratic
polynomial over GF(16), evaluate its cube at seven points, and interpolate the
degree-six result. Each GF(16) cube uses the existing four-AND cover; all basis
changes and interpolation are binary linear. Exhaustive algebra checks cover
all 4096 elements.

`check_gold_kernel.py` verifies all inputs and arbitrary target words under two
measurement seeds: 8192 shots, 28 Toffolis, 27 qubits including operands.
`check_controlled_gold.py` additionally verifies both coherent control values:
16384 shots, 40 Toffolis, 41 qubits. Controlling the source before cubing avoids
the extra cost of controlling each linear correction separately.

`gold_arithmetic.py` recognizes Frobenius-related products, prices the existing
terms, and uses the specialization only when cheaper. The first naive
controlled implementation increased the core to 800940 Toffolis; it is rejected.
The corrected source-masking version measures 796796 core Toffolis,
2072 below the selected DP core, at unchanged 1004-qubit width. Its complete build, default and independent native verification, and archive
audit all passed. `run_gold_candidate.sh` finished successfully. The selected
complete archive is frozen in `checkpoints/9.838916402692796`.

## Other independent research evidence

- `coefficient_rank.py`: the 47 coefficient bits plus constant and nine address
  bits have full binary linear rank 57. Simple exact affine coefficient-bank
  compression provides no saving for this table.
- `rank_trim_screen.py`: ten exact fixed-word chemistry rank/gauge screens
  achieve errors 19.55–19.91 mHa at C349–340, but effective normalization rises
  to roughly 57.5–57.9, versus 48.5673. These are numerical screens, not accepted
  submissions. Small rank savings do not compensate for that increase.
- `blocked_cf_model.py`: the blocked C-fraction tail identity passes 160 generic
  field/polynomial tests, including N70 blocks 4 and 8. Zero-pivot charts remain
  outside its scope. `profile_blocked_tail.py` counts the actual proposed tail
  arithmetic but deliberately excludes head extraction, transfer-matrix
  construction, determinant preparation, rational synthesis, reverse, and scan.
  Forward tail-only totals for block sizes 2, 4, 6, 8, 10 are respectively
  594830, 380767, 288094, 229222, 188055 Toffolis. These are not full estimates
  or scores; head cost and chart handling remain significant.
- `block_dictionary_screen.py` tests a different chemistry shape route: split
  the current exact-word dictionary into two blocks, then solve the unrestricted
  block-covariance span fit analytically. All six fixed partitions fail the unrestricted span screen: errors are
  281–1274 mHa, with normal-equation residuals below 5e-14. This excludes
  those fixed partitions even before PSD/rank conditions, and excludes no
  repaired or differently constructed dictionary.
- ABC truth-table synthesis of the rounded 9-bit square took 531 ANDs. This
  particular truth-table implementation is rejected; arithmetic squaring and
  other representations remain open.

## Reproduction

From this directory:

```bash
source env.sh
"$EXP_PYTHON" check_dp_arithmetic.py
"$EXP_PYTHON" build_dp_decoder.py
flock /tmp/r6-fullverify.lock "$EXP_PYTHON" \
  ../../../from_host/verifiers/dqi/verify_submission_dqi.py \
  candidates/gap_gold_n70_m4095.kmb --max-workers 3
"$EXP_PYTHON" gold_interpolation.py
"$EXP_PYTHON" check_gold_kernel.py
"$EXP_PYTHON" check_controlled_gold.py
```

Packaging uses the unchanged round-five integration packager, with local output
and local DQI evidence, chem-census `structured` evidence, and explicit reuse
of all three fully verified task reports. `results/package_dp.log` records its
completed audit. Re-run `bank_milestones.py` after each accepted package.

All modifications and generated artifacts are under `initial_experiments/6/agent14`.
Earlier experiment sources and host code are read-only dependencies. Their
Apache-2.0 attribution is retained in `initial_experiments/3/import_published.py`,
its pinned vendor tree, the round-five reports, and the bundled host licenses.

## Quadratic-extension Schur architecture (components only)

`quadratic_schur_model.py` moves interpolation nodes outside GF(4096),
retaining all magnitudes. A reversible affine frame reduces each remaining-node
update to an in-place GF(2^24) inverse and binary-linear transforms.
`extension_inverse.py` implements that inverse in 300 Toffolis, 97 total qubits;
32,768 native shots pass, including both axes, random extension elements, zero,
repeated composition, restored scratch, and phase. Private classical outcomes
are reset at component boundaries; this is necessary for repeated use.
`schur_pivot.py` passes forced zero-pivot tests in both directions.

`schur_extract.py --n 70` passes 256 native shots in each direction, split
between arbitrary data and actual syndromes of all weights, plus both
compositions. Full extraction costs 194,005 Toffolis in each direction and
925 qubits. This excludes polynomial synthesis, phase scan, and chart repair,
and is not a complete decoder or a submission score.

The independent C++ distribution census (`schur_census.cpp`, seed 142041)
with adjacent swaps recovered locations and magnitudes in 99,995/100,000
cases. Five extension poles remain unsupported; 786 swaps were used.
Classical swaps are not yet a quantum circuit.

## Completed narrower inverse and chemistry trials

The curvature-7 inverse fit passes all 3,133,452 integer inputs with 324 leaves.
`build_narrow_inverse.py` compiled nine schedules; the best Q88/T1863 objective
185,533 exceeds the selected 182,956.5. It was not promoted. Narrower slope or
constant precision failed their numerical constraints.

`robust_fixed_fit.py` completed 12,000 steps from the B511 dictionary. Final
continuous error was 48.949 mHa; rounded error 175.898 mHa, and rank-refit
132.582 mHa. `quantized_fixed_fit.py` then trained with exact rounded forward
words and a straight-through heuristic derivative for 16,000 steps. Rounded
error was 190.173 mHa and refit 166.686 mHa. Both fail the error limit and are
rejected. These failures concern the tested schedules and do not establish
impossibility for the shape.

## Complete quadratic Schur decoder and coherent repair

The earlier component-only limitation has now been resolved into a complete
native-verified phase oracle. `schur_synthesize.py` composes degree-by-degree
proper/monic rational synthesis, overwriting the consumed 840-bit pivot bank.
`extension_division.py` repays measured pivots through R(alpha)=a+b*alpha;
its create and parity kernels pass arbitrary extension operands, including zero.
The base two-word synthesis is a bijection even on a zero divisor.

`spooky_extension_inverse.py` delays the forward old-input HMR phase until
reverse inversion reconstructs those controls. Its forward/reverse costs are
192/300 T, with 16,384 native pair shots including both subfield axes and zero.
`schur_extract_spooky.py` passes full-size arbitrary-input tests, predicting
its pending forward phase exactly and clearing it in the paired inverse.

A complete un-repaired decoder costs Q949/T885278, passes the full verifier,
but has 756 value errors and 1 phase error. Repetition makes its cost
1641213135806.1492 worse than the earlier selection. It was not promoted.

`schur_extract_repair.py` coherently swaps adjacent node samples when the
current actual b vanishes and the next does not. Consecutive swaps are suppressed.
All 34 node descriptors stay quantum. Controlled constant field scaling needs
only 12 T: create the selected linear output, HMR the old word, and repay
against the inverse selected map. Exhaustive source/control checks pass.
`schur_dynamic_fold.py` uses a structured factorization of the polynomial D
map rather than controlling a large dense binary map. Both norm and quotient
phase are adapted to the coherently selected irreducible quadratic.

Full extraction tests include forced swaps, all true error weights, arbitrary
syndromes, exact pending ghost phase, frame/node descriptors, inverse, and clean
paired restoration. Full synthesis tests include 1,028 set node controls across
256 cases, all coefficient values, both directions/compositions, scratch/phase.
Selected forward/reverse extraction costs 175440/239700 T, and synthesis
137352/141432 T. `schur_scan.py` supplies the remaining 327600 T at 80 T/point.

The native scan uses monomial zero padding for completed tails, includes the
implicit monic leading coefficient, and divides the Forney residue by a*P(a).
A missing newline in the development test header initially dropped the first
basis-change gate; exhaustive testing caught it. The corrected test, small
complete oracle, and both full verifications pass. No verifier was modified.

Reproduction: build/test repaired extraction and synthesis first, then
`build_schur_decoder.py --n 70 --m 4095 --repair`, followed by the unchanged
full DQI verifier. Complete build and both native reports are terminal.

## New fused inverse arithmetic (not yet a full candidate)

`rounded_square.py` now constructs h=round(u^2/64) for every u=0..256,
including the endpoint guard. Split u=8v+r and use
h=v^2+floor((v*r+k(r))/4), k(r)=floor((r^2+32)/16).
A shared two-bank ANF evaluates the five-bit square; the low two carry bits
are measured and repaid against a five-bit source function. Forward/inverse
cost 77/83 T at Q37 including operands. All 257 magnitudes, two seeds, both
fresh directions/compositions, phase, and scratch pass (4096 shots).
This replaces the earlier losing 531-AND truth-table implementation, not an
unbuilt optimistic estimate.

`rounded_coordinate.py` prepares abs(floor((z+16)/32)) from centered 14-bit
mantissas with one's-complement folding, including its rounding carry guard.
It then calls the square and restores the mantissa. Costs are 93/99 T and
Q51 including 14 input and 11 output qubits. All 16,384 mantissas pass in each
direction plus both compositions (32,768 shots), with exact phase and clean
scratch. The source-independent extra peak is 26 qubits, including the
9-bit rounded coordinate; clean phase headroom can potentially be borrowed
in the complete circuit, but such a saving has not been measured yet.

`fused_round_screen.py` proves once-rounded full linear/curvature products
pass even without an extra fractional output bit, using the existing table.
That arithmetic still needs efficient truncated products. More directly usable,
`fused_booth_screen.py` models the actual streamed Booth/binary row rounding,
then refits only an integer c0 bias per existing prefix leaf. With q=5 and one
extra fractional output bit, both Booth and binary curvature rows pass all
3,133,452 inputs with 269 leaves. The Booth table has worst ratio
0.8454602023035518 and minimum phase margin 119742674 units. Its local table is
`parameters/fused_q5_e1_booth.npy`. This is a numerical model, not a full
native circuit or score. The full circuit must price all coefficient replay,
live registers, carries, shifts, rotations, and cleanup.

`fused_leaf_repair.py` tries +/-5 integer slope/curvature changes within each
old leaf to eliminate the extra fractional bit. Some leaves remain infeasible;
those tables are rejected. No model FAIL is a candidate. These are finite
family failures, not a lower bound on all fused approximations.

## Further fixed-dictionary chemistry screen

`overlapping_dictionary_screen.py` admits shared atoms in exact-word R2
blocks, beyond the earlier disjoint partitions. Twelve dictionaries use B384,
400,416,448, with 49,81,113,177 shared rows chosen by hub, weight, or fixed random
criteria and covariance spectral partitions. Every unrestricted, possibly
indefinite/arbitrary-rank block-covariance optimum fails: error 189--868 mHa.
The Moore-Penrose shared-intersection chart has normal-equation residuals below
3.5e-14. This excludes these fixed dictionaries only, before PSD/rank/physical
conditions. It does not exclude repaired geometries or new factorizations.

## Architectural notes to continue

1. Direct synthesis still uses four field products per coefficient through
   S_i D_i S_i^-1. Reparameterizing consecutive S_i^-1 S_(i+1) into affine
   factors can halve those products, but in-place source erasure and endpoint
   normalization must be built, not assumed.
2. At most one pivot swap could be recorded as a six-bit index instead of
   34 bits, with rare additional approximation failures. Measure that law
   independently before accepting a resource/error tradeoff.
3. The fixed-node scalar chart has now passed independent classical algebra
   and native extraction tests. See the new section below; it supersedes the
   earlier orientation of the degree-two transfer written during derivation.

## Fixed-node scalar chart: latest uncompleted architecture

`schur_gap_model.py` recovers all locations and magnitudes in 216 N6, 216 N14,
and 720 N70 cases, including 340 deliberately forced initial scalar pivots at
N70. It stops completed zero tails classically and is not a full circuit.
`schur_gap_census.cpp` independently samples the true DQI weight law: seed
142060 passes 99,999/100,000 cases, with one extension pole and no unsupported
scalar charts, false completions, or location/magnitude errors on supported
cases. There are 725 scalar chart entries. The independent fixed seed 142061 passes 99,998/100,000: one extension pole
and one non-completion zero value, with no supported-case location/magnitude
errors. Both runs are terminal. These are model censuses, not native full scores.

`schur_gap_pivot.py` implements the changed affine frame. Its 12,288 native
shots include arbitrary pivots and the current/previous/prior zero patterns,
three boundary charts, both directions/compositions, restored scratch/phase.
`schur_extract_gap.py` composes it with the delayed-repair inversions and fixed
node updates. Full N70 extraction passes 256 cases per direction, including
arbitrary data, completed tails, actual syndromes of all weights, and forced
scalar pivots, plus exact forward ghost phase and paired cleanup. It costs
132123 T forward and 196383 T reverse at Q921/925, with no separate node flags.
This saves 86,634 T versus the selected repaired extraction alone. **Synthesis
for this representation remains unbuilt**, so this is not a new full score.

Correct algebra: write q_i=x^2+x+k_i. If R(alpha_i)=a is scalar and
z=(k_i+k_j)/(R(alpha_j)+a)=u+h*alpha_j, define v=u+h*x. The synthesis matrix is

    M_gap = [[q_i+a*v, a*(k_i+k_j)+a^2*v], [v,q_j+a*v]].

Its determinant is q_i*q_j, and its inverse tail map is

    R_tail = a + q_j / (v + q_i/(R+a)).

The raw reciprocal updates therefore retain the **original fixed node
constants**. If S=[[a,1],[1,0]], S^-1 M_gap S=[[q_j,v],[0,q_i]]. The matrix with
the q_i/q_j diagonal reversed is the forward decoder matrix, not synthesis.
The independent model checks this orientation against all true magnitudes.

Equivalent serial synthesis factors are:

    normal: [[x+a/b, (a^2+a*b+k_i*b^2)/b], [1/b,x+a/b+1]]
    gap_first: [[a,q_i],[1,0]]
    gap_second: [[u+h*x,q_j+a*(u+h*x)],[1,a]]

The second factor uses the preceding scalar a. After gap_second acts on a
proper degree-d tail, the intermediate numerator has degree d+2, denominator
degree d with leading coefficient a; gap_first makes the pair proper/monic
again. The mixed coefficient layout is the next concrete implementation task.
For k=d+1 consumed stages, the pending state can store N0..Nk and B0..B(k-2)
in exactly 2k field words: N(k+1)=1 and B(k-1)=a are implicit. Use the last B
slot for the extra N coefficient. Only a temporary field may be needed when
restoring the ordinary coefficient layout. Completed a=b=0 tails can be padded
by x as in the proven proper synthesis; classify them before the scalar chart.
Do not assume a degree-mixed reversible layout or its resource cost until built.

A possible common factorization is M=S_left D_i S_right with
D_i=[[x,k_i],[1,x+1]] and S_left=[[b_arm,a+g1*x],[0,1]], where b_arm=b or 1,
and g1 indicates a non-completion first scalar pivot. The right factor is
[[1/b,a/b],[0,1]] normally, [[1,x+1],[0,1]] for gap_first, and
[[1,a_previous+(x+1)/b],[0,1/b]] for gap_second. For the latter, scale B by 1/b
first, after which the x-shear coefficient is the Boolean chart flag; this
avoids a fresh quantum field multiplication per coefficient. Conditional
field routing and implicit leading coefficients still need native tests.

## Latest inverse primitive refinement

The shared low/high ANF banks reduce the square to 77 T forward and 83 T inverse
at unchanged Q37. The low-carry phase uses only two ANDs plus conditioned CZs.
The full mantissa-to-h primitive is now 93/99 T, Q51; all earlier exhaustive
native tests pass again. `rounded_square.py` and `rounded_coordinate.py` contain
the selected component implementation, not the earlier 97/113 T drafts.

The extra=0 +/-5 per-leaf slope/curvature refit leaves unsupported intervals
and is rejected. The existing extra=1 numerical tables remain passing and
retain all 269 lookup leaves. No full fused inverse circuit has yet been built.

## Native two-node scalar synthesis component

`schur_gap_fold.py` builds the degree-two transfer directly in the existing
coefficient bank plus four fresh words, with its leading coefficient implicit.
It applies S^-1, the triangular q_j/v/q_i map, then S. Old head words (a,0,u,h)
are HMR-erased sequentially: recover (u,h) from
(k_i+k_j)*B(alpha_j)/(N(alpha_j)+a*B(alpha_j)) while a remains live, then a
from R(alpha_i). Fresh inverse creation uses the same two exact identities.
The top implicit monic coefficient is initialized in the fresh bank before
applying the invertible truncated polynomial map. Completion zero heads pass.

Native checks at tail degrees 0, 2 and 33 all pass 512 shots per direction,
including arbitrary head fields, monic proper tails, zero completion cases,
full coefficient-model equality, both fresh directions/compositions, scratch
and phase. At d33 the component costs Q961/T6390 forward, Q969/T6630 inverse.
It is not a coherent dispatcher between normal, scalar, and completion shapes.
Mixed-degree synthesis across a scalar pair straddling a fixed block boundary
is still required. **Do not replace the selected complete decoder with it.**

A useful allocator detail was resolved: Kickmix public registers created by
create_quantum_register are reserved and cannot be freed. To recycle erased
input head words inside this component, allocate their qubits normally and
attach register declarations in the final text. Both native directions and
composition tests validate that reused public head slots finish with the
required value. Applying this technique to the old proper/dynamic fold may
reduce its forward peak; its reverse creation peak then becomes relevant.
Any such rewrite needs native component and full archive verification.

## Concrete affine-chain route (unbuilt)

For nonzero pivots, S_i=[[b_i,a_i],[0,1]] and
M_i=S_i D_i S_i^-1. Reparameterize the original heads backward into outer
(a0,b0) and T_i=S_i^-1 S_(i+1)=[[t_i,u_i],[0,1]], where
t_i=b_(i+1)/b_i and u_i=(a_(i+1)+a_i)/b_i. Each T only scales N and adds u_i*B,
using two field products per coefficient. At a fixed quadratic node after D,
the ratio is alpha, so after T it is u_i+t_i*alpha; the tested extension
quotient phase can erase those source parameters. This makes a native T-step
component a useful next prototype. The intermediate numerator and denominator
have the same degree: store one extra numerator leading field. The final outer
S0 makes that leading numerator zero, allowing it to be freed.

Completion must be included explicitly. Arm each original zero b to one and
keep the count of the non-completed prefix in six qubits so inverse
reparameterization restores zeros. A completed D should multiply by x, not
introduce k/(x+1). On an all-zero completed rational tail, the ordinary D map
needs only its two top constant contributions removed, controlled by i>=degree.
For T after such a D, the high component of the quotient is zero. Restore t=1
for a completed T, or t=1/b_last_nonzero for the boundary T, and repay exactly
that omitted HMR phase. An extra prefix-product field, possibly reusing the
old extraction frame after its HMR, can supply the boundary factor. The prior
source/leading fields and node flags affect width; full cleanup and count work
must be included. This remains a design, not a claimed 100k saving or score.

## Affine-chain implementation and scratch contraction

`schur_affine_parameters.py`, `schur_affine_dynamic_step.py`, and
`schur_affine_synthesize.py` implement the earlier affine-chain design.
Both complete synthesis directions and their composition pass native N6/N70
checks on supported charts, including all prefix lengths and variable node
descriptors. The complete oracle passes default and independent full runs.
The public frame is measured in X, its phase is held explicitly, and its
inverse recreates the original frame and repays that phase.

`schur_extract_index.py` records at most one adjacent swap in a six-qubit
index. A later scalar pivot is unsupported and contributes to full verifier
errors. Local node controls are decoded from the index and cleaned after use.
Both extraction directions, forward ghost phase, and composition pass native
N6/N70 tests with forced swaps, arbitrary inputs, and actual syndromes.

`schur_dynamic_division_compact.py` forms the quotient numerator directly in
its output, then scales it in place. The inverse creation uses fewer scratch
fields at unchanged T. `schur_dynamic_division_contract.py` streams runtime
classical masks into one contraction field using the native bilinear MAC,
instead of materializing two numerator fields plus a contraction field.
The phase cost stays unchanged, and native d3/d34 tests and full synthesis
checks pass. Full DQI width falls from985 to955 at unchanged951995T.

`lifetime_color.py` found only one saved qubit in the large forward affine
step and none in reverse. That allocator was not integrated. The alternative
`schur_remainder_step.py` and dynamic counterpart pass native component
checks, but coherent node control raises the dynamic variant's T too much.
They are not selected decoder components.

## Completion frames and seedless synthesis

`schur_affine_parameters_filled.py` fills zero suffix heads with the last
active affine frame. Completed T factors become identities, and completed D
steps multiply the constant rational tail by x. The old boundary reciprocal
and its flag are removed. The old measured frame's twelve wires can hold
scratch. Full synthesis N6/N70 checks pass at Q940/942, T98492/108495.
The complete default oracle with the original separate scan scratch (Q947)
passes102400shots, 8 value errors and1 phase error. The Q942 scan-reuse build
has not received full evidence: its attempted run failed at import because
that shell omitted env.sh. Do not reuse the empty JSON as verification.

`schur_affine_parameters_canonical.py` and
`schur_affine_dynamic_step_canonical.py` remove the degree counter. Zero-head
predicates are measured in X and their private outcomes retained until the
inverse. The inverse recognizes copied heads by equality; each step recognizes
completion by an identity factor and by a zero high quotient on reverse.
This is an approximate chart: equal nonzero active heads and accidental active
identity factors are unsupported. The full verifier charges their incidence.
Native N6/N70 supported-chart checks and the small complete oracle pass.
The Q936/T954283 complete default passes102400shots with4 value errors and
0 phase errors, cost1011164580515.3181. It is not yet independently checked.

`schur_affine_parameters_seedless.py` puts the endpoint ratio in the original
outer a0 slot through an invertible affine Horner map. The custom
`schur_affine_outer_seedless.py` consumes b0 and the old numerator leading
ratio, erasing their phases through the two components of N*conjugate(B).
It removes the extra endpoint field. Native N6/N70 synthesis checks pass at
Q924, T104430/115159, including pending completion/frame ghosts and paired
cleanup. It inherits the canonical chart's explicitly unsupported cases.

`spooky_extension_inverse_compact.py` allocates its output after the norm
inverse scaffold has been released. Reverse HMR source slots can host phase
scratch. The full inverse pair passes16384native inputs, axes, zeros, random
values, held-ghost phase and cleanup at Q75 in both directions, unchanged
192/300T. `schur_extract_index_compact.py` passes N6/N70 extraction checks at
Q912 in both directions, unchanged176460/240720T.

`build_schur_affine_seedless.py` composes these components and reuses only
known-zero frame, source, and dropped-leading wires for scan scratch. Its
native N6/M63 oracle passes4096shots with no errors. The default full report
and ongoing independent report have the prefix `schur_affine_seedless_n70`.
Only promote after digest-bound packaging and audit. Goal14 remains active.

## Milestone11 acceptance and independent canonical census

The seedless default and independent full runs both pass. The exact archive
was packaged, re-extracted, digest-checked, audited, and frozen at threshold11.
`schur_canonical_census.cpp` independently checks two fixed100000-shot
host-law streams (seeds142078,142079): respectively99995/99992 supported
passes,4/4 pivot exceptions after the one-swap limit,1/4 extension poles,
0/0 false completions,0/0 supported locator or magnitude errors. The streams
contain1/0 coinciding-head canonical ambiguities, reported separately. These
are algebra censuses, not native score samples. All pipelines are terminal.
Goal14 remains active, and thresholds12,13,14 remain outstanding.

## Post-milestone13 width and converter changes

`inverse_ratio_mask_phase.py` passes16384 native cases at Q47/T69, using a
classically factored scale on the divisor and the tower inverse-dot kernel.
`derivative_inverse_deferred.py` keeps companion HMR masks across the phase-only
middle scan. N35 native256 including forced top gaps PASS, Q883,
T80475forward/82726inverse, saving12Q and3541T versus the raw converter pair.
`cubic_mask_phase_compact.py` fuses GF16 inverse output linear forms directly
into the trace CZ pairing. All57344 native cases PASS at Q27/T30 (four fewer
scratch qubits than joint kernel). Not yet integrated into a full submission.

`raw_cf_endpoint_seedless.py` transforms the first source head directly into
the projective terminal jet, no copied tail fields. NativeN70 PASS512, Q867,
T6405eachdirection. `raw_cf_outer_seedless.py` reconstructs the missing outer
head from the leading pair, contracts the two known-zero top coefficients,
and repays the outer head ghost on inverse. NativeN70 PASS256, Q901/902,
T1777/2328. Full synthesis seedlessN70 PASS256, Q902/903, T39905/59190.
Complete small N6/M63 PASS8192 with no failures. Full seedless Q903/T554266
(default6value2phase, cost552738368374.0004; independent1value0phase,
cost540860096491.75934, saltagent14-raw-cf-seedless-independent-a) was promoted after
all eight fixed default failures matched the old raw circuit in two regression runs.
No hash/salt tuning; existing milestone13 archive remains selected.

Chemistry incumbent B671 ADMM terminal: best integer refit remains worse
than initial23.943mHa; endpoint49.316mHa, invalid. A different integer-only
1000-round descent now runs from the exact incumbent B671 dictionary. Its
first185 rounds reduce23.943->21.282mHa, still invalid. C340 forty L1 Jacobi
sweeps reduce effective normalization only57.9->56.125, insufficient to offset
small rank savings; no such chemistry fit is promoted.


The full-angle projected-Hessian integer method (`integer_newton_fit.py`)
accepts1652 row moves over150 steps from the B671 greedy endpoint, reducing
21.127105->21.105792mHa. Still invalid; lambda_eff56.30959. This is scoped
negative evidence for the correlated QP proposal family at this seed, not a
shape-wide bound. `rank_one_blocks_screen.py` direct physical-factor splitting
intoR350/C1 gives float errors4.771Ha atB7 and2.985Ha atB15; the simple spectral
truncation initializer is not competitive. Exact B76 reproduces the float
incumbent error before rounding, confirming the split implementation.

`quartic_mul_seed.py` builds a distinct GF8 quartic 42-AND multiplication seed
(contrast upstream GF4 CRT). Algebra28672cases PASS. The corresponding local
quadratic-cover walk retains42; common-operand two-GF8-product search retains12.
These are search outcomes, not optimality proofs. `schur_division_lowtemp.py`
recycles output slots and erases the copied numerator against untouched inputs.
Native16384 PASS including both axes and zero, Q94/T712 (74 public qubits,20
private), trading extra gates for6 fewer scratch qubits than cubic quotient.
No complete candidate uses that trade yet.


`raw_cf_synthesize_lowtemp.py` additionally clears the padding flag using the
padded denominator constant term (all actual roots are nonzero). The outer
head computation lends a zero head field to normalization and then measures
and recreates the other head to restore its borrowed denominator coefficient.
Its inverse erases copied heads directly against their unchanged sources.
N70 outer256PASS Q895/896,T1902/2440; full synthesis256PASS Q895/896,
T40041/63155. Latest small complete8192PASS, zero failures. CompleteN70
Q896/T557645, SHA7cc5863d695ab940fd9f4190367b4c3b09fa8ac718ce7a848d56a35be708ada2,
full default pending. No score improvement claimed before native evidence.
The previous13.419267220752062 archive remains selected.

A first coefficient/master-angle STE chemistry fit (with whitened factor
weights) currently degrades the B671 integer tensor. Best checkpoints retain
the original state. A separate projected-STE implementation is prepared: its
integer forward objective analytically minimizes the PSD rank-C coefficient
fit, removing stale-weight effects from angle proposals. It is not run yet.
Both use exact integer native tensor checks for acceptance, regardless of the
heuristic backward derivative. New inverse address-embedding coefficient fit
with three prescribed C0 bits misses the bound (6.08967e-8 vs5.96046e-8);
no such coefficient table or circuit is promoted.


Lowtemp full default PASS102400:6value0phase, Q896,T557645max,
557637.9150488281mean, cost544628838530.43414. Fixed independent salt is
agent14-raw-cf-lowtemp-independent-a, currently running. Common8seedlessfailure
cases reproduceexactly under2measurementseeds, supporting unchanged arithmetic
semantics. Promotion awaits the independent report and archive audit.


Lowtemp full independent PASS102400,7value0phase, cost546103877154.05615;
package+ZIPreopening+digestbinding auditPASS13.485543697862889. Frozencheckpoint
asabove; milestone13 unchanged13.202469964520382. Only14 remains.
STE6000stepsandprojectedSTE6000stepsterminal: neitherimprovesinitial21.105792mHa;
endrefits116.387mHa and102.230mHa respectively. Bothbestfilesretaininitialstate.
A fresh native-metric LLLproposalfit nowruns fromthe21.105792mHa B671state,
using fullprojected75-angleHessians, trust64, and exactintegerrefitacceptance.
Prior4latticecodeisreadonly; its localbarriercontrolreproducesexhaustiveoptimum.
Upstreamandigu/uncompute HEAD was checkedread-only andremainsthesame frozen
64b5a13aa249f1360846657392e6c66817b26d0f; no newexternalcircuitwas imported.
