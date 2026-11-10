> Current handoff: [REPORT.md](REPORT.md), written 2026-10-25 after the capacity interruption. Packaged best is **12.046515352810363**, with milestones **10–12** saved. A newer recovered-gap decoder passed both full native checks and awaits packaging; its combined verified-cost projection is **12.11607972269253**. Status statements below are historical and are superseded by REPORT.md. No experiment jobs were running at the handoff check.

# Active objective

Build a fully verified combined submission scoring at least20, saving archives
at every integer threshold10–20. Goal remains ACTIVE. Do not stop at milestone10,
reduce scope, change verifiers, or search salts/hashes. No subagents are allowed
by the current developer instruction. All new writes stay in this agent20 tree;
prior experiments and host verifiers are read-only.

## Accepted repaired-H decoder — latest work

Both fixed full native checks and the original combined metric PASS. The new
accepted `best/submission.zip` scores **10.128965274112238** by default and
**10.144446072819372** with the predeclared independent changed-task salts.
ZIP SHA256 df518f8d20a0b4d4ea72700280f689bc17d0685471e0d1ddb2cc998b1c60a52e.
`results/accepted_H_improvement.json` records promotion; original metric audit
is `results/audit_best_H_metric.log`. Previous best is preserved at
`verified_before_H`; the first score-10 archive remains immutable.

The repaired-H endpoint is Q1004/static1123942, versus old1162746.
`candidates/repaired_H_dqi.kmb` SHA256
ac2d7a963c1b3e2f58f39a08000838296f4a6fd18d468cd35567351b211d41ed.
Default102400 shots: PASS, zero value/one phase error, cost1355783099065.6162.
Independent102400 shots: PASS, zero value/phase errors, cost1349585429383.3171.
Independent salt `agent20-repaired-H-dqi-independent-20261025` was fixed before
runs. Reports `results/dqi_H_{default,independent}_verify.json`; timing and
hash evidence `results/dqi_H_verification_runs.json`. No verifier is running.

Source/byte manifests: `parameters/repaired_H_pre_native.json`,
`results/repaired_H_build.json`. Sources `H_decoder.py`, `repaired_sylvester.py`,
`fast_sylvester.py` built and verified these bytes. Small complete N14/M127
oracle PASS4096 shots, zero value/phase errors; component H degree7/17/35
PASS1024 shots each including forced first gaps/all padding degrees and
arbitrary paired inverse; native-field law model50000 vectors has zero errors
with170 repaired gaps. Full core Q1004/T925948; endpoint46*4095+9624.
Reports `results/repaired_H_small_verify.json`,
`results/repaired_sylvester_H_checks.json`, `results/H_decoder_model.json`.
Goal20 is still ACTIVE. Milestones11–20 remain missing.

Latest numerical chemistry screen COMPLETED, no GPU job live:
`chem_residual_grow.py` grows R4/B224 exact rounded words toB257 using missing
generator-gradient singular vectors, quantized matrix eigenvector atoms and
the accepted dictionary pool, with exact PSD block refits. Error153.94→84.60
mHa, Lambda1055.29,lambda_eff94.11: fails accuracy and normalization. Nine
finite growth batches, no restart justified. Report
`results/residual_grow_R4B257.json`, checkpoints parameters sameprefix.

## Authoritative accepted submission

- `best/submission.zip`: SHA256
  df518f8d20a0b4d4ea72700280f689bc17d0685471e0d1ddb2cc998b1c60a52e.
- Full local default combined10.128965274112238; independent changed-task salts
  give10.144446072819372. No Kaggle upload has been made.
- Inverse `candidates/rows_odd_best.kmb`, SHA256
  ddc952552caef19c5645ee4dcfc6c1b244d2353f6636d9b40118a458dc155105,
  Q89/static1709/34 rotations, cost173447.66648015822. Both200540928-shot
  native checks PASS (independent cost173447.6937108191).
- DQI `candidates/repaired_H_dqi.kmb`, SHA256
  ac2d7a963c1b3e2f58f39a08000838296f4a6fd18d468cd35567351b211d41ed,
  Q1004/static1123942; full reports and costs above.
- Chemistry unchanged round5 `chem-census/candidates/structured`, full50880
  shots, cost696942846.9115301.
- `best/results/archive_audit.json` and `metric_preflight_and_score.json` bind
  every submitted byte to native evidence and use the original metric.
- `milestones/10/submission.zip` preserves first-crossed10.030435786719812,
  SHA25659bf461be714e0a4786fb0cc825d8b105622fec35ce193acd44ab1298004c2a0;
  previous best directory is `verified_before_H`.
- `verified_inverse` preserves intermediate9.98716898063121 submission.
- Milestones11–20 remain missing.

`package_best.py` uses original round5 packager with all destinations local and
saves future crossed thresholds10–20. `audit_best_metric.py [directory]` uses
the original readers and score; invoke with chemistry Python for pandas.
Do not rebuild accepted bytes for a tiny projected gain.

## Environments

`source env.sh`; EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python (3.13 native
kickmix). Existing imports include round4/codex, round3, round2, host. All local
source copies are here; upstream andigu/uncompute Apache2 rev64b5... is read-only.
`source chem_env.sh`; EXP_CHEM_PYTHON=/mnt/data01/miniconda3/envs/torch2.13/bin/python.
CUDA RTX4090/24GB. GPU jobs hold `/tmp/r6-gpu.lock`; full native verifiers hold
`/tmp/r6-fullverify.lock`. Cap own CPU modestly (8 shared CPUs). Cache/temp
writes local. python-flint0.9 installed locally under `dependencies` for an
unsuccessful symbolic factor screen; prepend that path only when needed.

## Completed finite chemistry pipeline — do not restart

Exec session57141 has finished; no GPU job remains from this pipeline.

`source chem_env.sh && flock /tmp/r6-gpu.lock "$EXP_CHEM_PYTHON" chem_alternating.py --init parameters/tangent_eff_B703_last.npz --cycles 5 --tag alternating_B703 > results/chem_alternating_B703.log 2>&1`

It alternates up to5 cycles of300 exact integer-word steps (normalization
penalty1e-4) and1000 preconditioned physical-weight steps, with Lambda<=330,
lambda_eff<=51, gap>=3.91504. The parent holds the lock; children do not relock.
Child logs:
`results/alternating_B703_<cycle>_<words|weights>.log`; histories same prefix;
checkpoints in parameters. Parent progress JSON is
`results/alternating_B703_progress.json`, written after each completed phase.
It preserves the minimum-error checkpoint at
`parameters/alternating_B703_best_error.npz` and records source/checkpoint hashes.
It stops early if numerical error2<.019995, reporting POTENTIAL_TWO_BODY_PASS;
this is NOT a native physical acceptance. Independently assess actual words,
one-body reconstruction, physical constraints and ANGLE circuits before any
submission selection. Do not edit its two child source files during this finite
pipeline without a clear reason (source hashes are recorded per phase).

Final best numerical state `parameters/alternating_B703_best_error.npz`,
SHA256 c3094a02b0591b33b20e1f1ee38f531595e8679ebc7c0fb6cedeaf2f91b4e7fd,
error20.144102913mHa, gap4.00263463, Lambda326.9123108,
lambda_eff50.99999990. It still FAILS20mHa. All five cycles completed.
Actual unchanged eight-bit words, tensor and physical checks are float64.

## Chemistry method and completed controls

`chem_fixed_weight_repair.py` ranks 1/2/4/8-angle row proposals with float32
exact-word gradients, then evaluates actual decoded dictionaries in float64.
It keeps the weight gauge fixed for word moves. Optional PSD-refit blending
plus normalization headroom first lowered B703 from20.552 to20.422mHa
(margin2), then20.414969mHa (margin1), but stalled.

`chem_tangent_refit.py` now uses actual residual gradient and a feature-Gram
preconditioner (ridge.001), projecting in that metric onto three linearized
physical halfspaces via an exact3-variable active-set QP: gap floor, Lambda cap,
and effective-normalization cap. Optional `--project-constraints` retracts
trial states to actual physical feasibility with float64 Newton corrections
and line search. Every accepted weight step actually lowers error and passes
all optimizer physical limits. These limits are search constraints, not changed
host acceptance rules. Offsets are projected toward native gap margin, but
sometimes the attainable minimum gap is higher.

Completed stages:
- tangent_B703: initial PSD-direction strict tangent control stopped step0.
- tangent_gradient_B703:500 steps,20.414969→20.314905mHa,
  Lambda319.1601,lambda50.8775.
- fixed_B703_after_tangent1:123 word steps,20.311776mHa,lambda50.8629.
- tangent_qp_B703:2000 steps, about20.295mHa,lambda50.917.
- tangent_precondition_B703:stopped117 at Lambda320 boundary,
 20.234433mHa,lambda50.9632.
- tangent_manifold_B703:stopped18,20.234316mHa,lambda50.9064.
- norm_words_B703:normalization penalty1e-4, stopped190 with5899 row moves,
 20.232904mHa,lambda50.7830 (some accepted word moves trade error for norm).
  `_best` is best objective; `_best_error` explicitly tracks minimum error.
- tangent_eff_B703:relaxes artificial raw cap320 to330 while keeping
  effective cap51, reaches20.1612016mHa. Live alternation starts here.
- Earlier projected QAT10000 steps and robust_banked31 8000 steps COMPLETED,
  no useful changed-word improvements. Do not restart them. QAT best feasible
  zero-changed-word error19.9528mHa; terminal changed-word95.79mHa. Robustness
  final rounded157.85mHa/refit51.98mHa. These are scoped negatives.

Exact rank-C capacity at B719 passes toC340 (19.9102mHa), but subset-aligned
physical gauges lambda~57.24 lose. PCA gauge~103 loses. `chem_rank_screen.py`
and `chem_rank_gauge.py` retain evidence; not a global constraint on new words.
Prior round5 R2/B351 continuous C189 fits19.98mHa, but rounded dictionaries
have140mHa span lower bounds. Existing block/PSD/discrete research is extensive;
read the reports before repeating it. The exact-word generator partition plus optimal Ky Fan column rotation was
implemented in `chem_exact_partition.py`: B351 best465mHa, B479 best237mHa.
Those finite screens fail; do not repeat them unchanged.

## Important accepted-source failure already corrected

Do NOT use `candidates/retained_inverse_dqi.kmb` without safe:50501 phase
errors/102400 despite0 value errors. Retaining norm intermediates through a
source/reciprocal swap is invalid. `safe_retained_inverse.py` excludes mutable
opening/gap-serial sites and reserves body scratch. Small N14/M1274096-shot
complete and both full accepted checks validate it. Accepted inverse-row change
is exact three-AND Karatsuba across adjacent radix-four products; phase/margin
unchanged, selected Q89 by resources. `row_search.py` reproduces it.

## DQI new research — not accepted full submissions

1. `cube_field.py`:28-AND GF4096 cube, seven GF16 evaluations, exact basis
   maps. Exhaustive4096 values/arbitrary targets/two seeds8192 native shots PASS.
   `cube_arithmetic.py` recognizes Frobenius-related products and controlled
   input masks. Generic serial Cfraction whole1376094/Q878 vs1450238 baseline;
   missing zero charts and not competitive.
2. `blocked_cfraction.py`:two-pair map with one long reciprocal. Independent
   native forward/inverse N6/14/34 PASS phase/scratch. Noncompact Q974/whole
   1201626. `compact_blocked_cfraction.py` reuses head:Q915/whole1223658;
   native tests pass. Unsupported zero charts differ from serial.
3. `jacobi_decoder.py`:in-register pair f=a/(1+bx+x^2g), normalizes all
   nonconstant coefficients by a, unit-series inverse preserves b. Synthesis
   is a shear/bank-swap/a-scale. Pair uses one scalar inverse instead of two.
   Block2 A=1+(b+d)x+(bd+c)x^2,D=1+bx,p=ac,
   P0=a*(b^4+c^2+c*(b+d)^2),P1=(bd+c)*a*(b^3+cd),
   g=A*T/(p+x^2 D T),T=A*R+P. Plain block2 Q938/whole1100320;
   `compact_jacobi.py` Q915/whole1104604. Native components N6/14/34 PASS.
   `check_jacobi_oracle.py` N14/M1274096 shots exact model contract PASS,
   but145 truth-phase errors on forced unsupported charts.
4. `jacobi_census.cpp` exact independent GF4096/host-law model:stream
   2026102701/20000 has121 wrong vectors for serial and blocked Jacobi.
   Expected random-mask phase incidence half; penalty erases resource gains.
   `general_jacobi.py` larger blocks generate short B/A and use a middle
   product for P, avoiding expanded symbolic P. Block3 Q987/whole1037982,
   block4 Q1035/whole1027914; native block3 N8/14/34 PASS model/phase/scratch.
   Still missing zero charts. Do not build full generic binaries for acceptance.
5. `lookahead_jacobi.py` repairs zero FIRST Jacobi pivot in block. Raw
   [0,u,v,w] encodes [0,u,h=w/u+(v/u)^2,e=v/u]. f=u*x/(1+ex+hx^2+x^3g).
   Its g=A*T/(u+x^3T),P0=u e^3,P1=u h(h+e^2). Generic and gap both return
   encoded head[a,b,h,e], keeping a unchanged so zero predicate is temporary.
   Native mixed-chart N6/14/34 pass with half fixtures forced leading-zero.
   Natural stream2026102801/20000 has75 wrong vectors vs161 serial on SAME
   stream: roughly half failures repaired; zero SECOND pivot still missing.
   Separate selected synthesis Q928/core1210382/whole1537982 loses.

## Shared synthesis breakthrough and current exact resources

`shared_jacobi_synthesis.py` replaces selected ordinary/gap polynomial folds
with a shared action. Conditional P/Q swap, a numerator shear, two scalar
scales, and rank-(n+1) polynomial product share both charts. New boundary fields
are computed before the tail. Old head/parameter words are measured early and
released, then their phases are repaid after reconstructing the original head
from low OUTPUT POLYNOMIAL coefficients. Do NOT replace this with Taylor
moments: common factors on armed charts carry extra information.

The local inverse uses N0=a,N1,N2,N3 and1,B1,B2,B3. Normal:
qj=Nj/armed(a),b=B1+q1,c=B2+b q1+q2,
d=q1+(B3+b q2+q3)/armed(c); encodedold=(a,b,c+bd,b+d).
Gap:u=N1,q1=N2/armed(u),q2=N3/armed(u),e=B1+q1,h=B2+e q1+q2.
Code saves q1 storage via b=B1+q1, computes c=B2+b B1+b^2+q2 and
encoded e=B1+(B3+b q2+q3)/armed(c), h=c+b e+b^2.

`compact_arming.py` retains one copied zero bit instead of the entire AND
ladder across scalar scaling; recomputes the ladder after restoring source.
Shared module installs it only into imported research Cfraction helper in its
process. All host and prior source files are untouched.
`variable_shift_inverse.py` coherently chooses r/(1+x^2r) vs r/(1+x^3r), by
one clean coefficient, conditional rotation, and controls only on extra output
writes. V2 aggregates selected output in a field word and controls each bit
once. Native N6/14/34 PASS both selector values, phase, scratch, involution.
V2 N34 Q436/T10447 (408 input bits plus selector).

Current full N70 EXACT RESOURCE PROFILE:
Q942, core984034, whole1311634 (80*4095 scan), extraction399545 each way,
synthesis100254 forward and84690 inverse. This still loses after remaining
zero-second-pivot error penalty; DO NOT build full for acceptance as-is.
Resource log `results/lookahead_shared_full_profile.json`.

Shared synthesis native N4/6/14/34 tests pass both directions, two seeds,
including half leading-zero charts and quarter zero-second-pivot common-factor
charts. Latest allocation scheduling was additionally checked in COMPLETE
N14/M127 oracle,4096 shots:0 value/scratch/model-phase/mask/twice errors,
0 truth phase errors on that stream (forced unsupported charts included).
Artifact `lookahead_shared/candidates/cfraction_n14_m127.kmb`, SHA256
 a567606ec6fa1d1ecd8a6763e2dc19ed6a52c14ae980b2b74ef11b95d397f59d.
Report `lookahead_shared/results/cfraction_n14_m127_lookahead-shared-oracle-a.json`.
Other earlier small binaries remain under `lookahead` or top-level candidates.

Potential next DQI work: a degree3/carry repair for a zero second pivot,
new low-width block-Schur/Hankel representation, or a phase-aware endpoint that
cuts the scan. Simply making current leading repair exact will not reach20;
need substantial algorithmic gate savings or other-task architectural gains.
Source-polynomial compiler evidence in `continuant_polynomials.py`/`field_dag.py`:
plain Cfraction b4 P needs253 field products; bounded LRU thrashes, tree b4
118kT. Not a viable expanded-head compiler. Short-series preparation is better.

## Newly completed rounded-square inverse baseline

`rounded_square.py` exact 56-AND integer network, bounded measurement-correct
cache. All512 totalized inputs/arbitrary targets/two seeds PASS. Cache32
cost98T/Q54 including20 publicbits; cache64 cost56T/Q78.
`fused_square_invsqrt.py` is NOW BUILT and exact native-model checked:
`candidates/fused_square_schedule_best.kmb`, SHA256
5748a098be62abfd23f06f2892f6e1fce39b0c066245b92836e6e2d459928218,
Q122/static2063/35 rotations, mean2051.96136. Two seeds,42334 distinct inputs,
84992 shots match independent `square_model.evaluate(fractional=1)` exactly,
all input/scratch restored. Full all3133452 numerical model passes ratio.8889,
margin4.2523e-12 turns. Retained-coefficient baseline LOSES (~282k vs173447),
so no full verifier run justified. This closes the earlier unbuilt uncertainty
for this schedule, not all possible fused-square designs. Reports
`fused_check_best/results/rows_exact_native.json`,
`results/fused_square_schedule_search.json`, `results/rounded_square_components.json`.

## New exact-word multiblock bootstrap screen — completed, fails

`chem_exact_partition.py` alternates joint OMP from accepted exact atoms and
the optimal balanced Ky Fan column rotation for two fixed projection spans.
This avoids quantizing continuous multiblock vectors. B351:3 seeds x4 cycles,
B479:2 seeds x3 cycles, each followed by exact rank175/block PSD refits.
Both fail badly (best errors around465mHa and237mHa, normalization poor).
Reports `results/exact_partition_B{351,479}.json`. No circuit build warranted.
Native target h1/G are dense, so there is no simple exact diagonal-sign
orbital symmetry in the supplied basis (h1 has all5776 entries nonzero).

## New reciprocal and polynomial-inverse primitives

`parity_series_inverse.py`: f=E(x²)+xO(x²), R=O/E,
G=E^-1/(1+yR²), H=RG. In-place O/=E; E=inverse(E);
E/=(1+yO²); O*=E. Native exact independent reciprocal, arbitrary ignored
constant word, phase/scratch/involution PASS N2/3/6/8/14/34, with and without
known-zero coefficient1. N70 Q843/T27762 vs old30158. Accepted DQI schedule
resource-only port saves ONLY616T (core790980/Q1004); no full rebuild.
Reports `results/parity_inverse_checks.json`, `results/parity_safe_profile.json`.

`polynomial_inverse_reversal.py`: first general in-place polynomial modular
inverse prototype by monic Jacobi-quotient reversal. Native degree1/2/3/7/17/35
PASS1024 shots each: arbitrary totalized states including zero pivots match
the independent model, phase/scratch/involution; modular-product identity holds
on all nonzero-pivot fixtures. Degree35 Q868/T114138 (840 publicbits).

`polynomial_inverse_bilinear.py`: better unnormalized quotient extraction.
For F degree n, G degree n-1, store q1=leadF/armed(leadG),
q0=(F[n-1]+q1 G[n-2])/armed(leadG). Tail R=F_low+q*G
uses ONE rank-(n+1) convolution rather than two scalar MACs per word.
Recursion (F,G)->(G,R) stores2 coefficients/step in consumed F slots.
After extraction F=c K(q), G=c K(tail q). Reverse q; multiply q_i by
c^(2*(-1)^i); invert c iff degree odd. Continuant reversal keeps F and turns
G into G^-1 mod F on generic charts. Armed zero pivots give an involution
on all states, but NOT modular inversion on degree-gap charts.
Native degree1/2/3/7/17/35 PASS1024 shots each, including25% forced zero leading
G. A later trim drops unused leading G from the convolution. Latest degree35
Q892/T73542, 852 publicbits (extra leading F word included).
Artifact `candidates/poly_inverse_bilinear_d35.kmb` (~223MiB due free CNOTs).
Latest SHA256 db6d7a798ecf7e65897f836f227665301fd6a580d4c430fcca45a6b2acd713f2.
Latest reports `results/poly_inverse_bilinear_trimmed_checks.json`; older
`poly_inverse_bilinear_checks.json` describes the untrimmed76398T version.

Important integration caveats: adding H adapter (~50040T) gives126438T per
direction BEFORE padding/degree repair. Ordinary zero-pivot failures (~d/4096)
may erase modest whole-decoder gate savings via R; no full candidate yet.
Common padding is also important: if B=P Btrue, H=P Htrue, H is NOT a unit
mod padded B. Simply making P squarefree does not fix the common gcd.

Variable-degree insight, NOW IMPLEMENTED AND NATIVE CHECKED: for B=x^m Btrue,H=x^m Htrue,
run the unnormalized extractor unconditionally for35 steps. After d=35-m
real quotients, the first fake quotient is (q0=0,q1=c), then all heads and
terminal constant are zero. Reverse/scale only the first d heads, using c
from fake-head q1 (or terminal constant if d=35). Inverse extraction should
return B unchanged and x^m Htrue^-1. No controlled long arithmetic is needed;
only prefix reversal and selecting c depend on existing completion degree.
Independent padded classical model passes1435 fixtures/all degrees0..35,
skipping5 interior-gap fixtures. `padded_polynomial_inverse.py` implements it:
rotate both quotient banks right by completion count m, then globally reverse,
which reverses the real prefix and moves the only fake nonzero coefficient
to the final head. Swap it to terminal c using [d<n], scale all head pairs
(suffix is zero), conditionally invert c on odd d, and route c back to
head[d].q1 with unary degree indicators. Explicit middle inverse reverses
every operation, including arbitrary unsupported suffix/flag states.
`check_padded_polynomial_inverse.py` latest native N7/17/35 PASS: forward all
degrees/padded regular pivots; paired inverse arbitrary polynomials AND
arbitrary non-unary flags; two seeds512shots each. N35 Q933/T78882 forward,
157764 paired. Public887bits =852coeffs+35flags. Reports
`results/padded_poly_inverse_trimmed_checks.json`. Both finite jobs completed.
At physical roots e=beta^(70-4m)*Gnew², so dynamic exponent correction is
still charged. A unary completion-dependent scalar is linear in flags; root
AND each of12 G-evaluation bits then CZ with flags implements that phase in
23T/point (11 root +12 pairing), rather than pretending a free11T scan.
Interior Euclidean degree gaps still need repair, perhaps a shared block3
covering quotient patterns111/12/21/3. Do not build a full acceptance candidate
until natural-law census and complete resources justify it.

## New most promising next H-scan lead: reachable monomial padding

`sylvester_padding_screen.py` independently mirrors the prior armed Sylvester
algorithm and checks H=sqrt((AB)') when A=x^m Atrue,B=x^m Btrue.
Contrary to the earlier generic-rank objection, it computes the correct H
on these RESTRICTED padded payloads:100 fixtures per degree0..35 have only a
few interior-gap failures, including correct cases with large common x gcds.
This suggests irreducible repadding is unnecessary for the46T H scan. The
accepted scan uses UNPADDED A,B with e=beta^70/(A B'). Therefore after H is
prepared in RAW padded banks, apply the existing monomial unpadding rotations
to H and B before the46T scan, then undo them. Do NOT use the same beta^70
mask on raw padded H: that misses beta^(2m). The previous-cofactor phase_shift
rotation commutes with unpadding and must happen before H preparation. Check
conventions against actual accepted compact-state frames. Hypothetical~1.08mT/Q1004 would
beat1.16m modestly ONLY IF interior Sylvester gaps are repaired.

NEW exact first-gap repair identity, classical model implemented:
In Sylvester recursion f=qg+r, if r[n-1]=0,r[n-2]!=0, replace r by x*r.
Write old c=gamma+x*c1 (gamma=c0). Recurse with c1 padded with high zero,
then add gamma*r to the low output. Source shift and payload shift are
conditional linear permutations. This is c*r=c1*(x*r)+gamma*r exactly.
Keep gamma in ONE12-qubit garbage slot across the H phase; reverse H later
restores c0 and clears it. Repair only first source degree-gap; retain17
one-bit repair tags (one per depth), so later gaps remain an O(n²/q²) error
budget. Source-common-x padding makes some payload c0 zero; no rank claim
is made on arbitrary payloads, but total forward/inverse must restore all.
Correction uses a12-bit masked copy of gamma at each level (repair flag),
then n-1 ordinary field products, ~6.4k extraT per direction over whole17
levels, plus conditional source/payload rotations ~several kT. Expected
H map~64k versus50k generic, still potentially saves~55k wholeDQI gates.

The V1 repair screen COMPLETED but exposed a padding error: gamma=c0 and
dividing c by x lowers its required common-x valuation. V2 FIX:
k=floor(m/2), gamma=c[k], replace c by c[:k]+c[k+1:]+[0], correction
gamma*x^k*r. This keeps k leading payload zeros. Screen
`results/sylvester_padding_repair_v2.json`:12/36000 failures. Classification
proves ALL12 have Btrue constant0 OR nonsquarefree Btrue, which physical
nonzero distinct error roots exclude. NO failures on valid squarefree,
nonzero-constant fixtures. See `results/sylvester_padding_repair_classify.json`.
Double interior gaps still need a larger natural-law census.

`fast_sylvester.py` already improves source division: a degree-one quotient
uses a single bilinear convolution. Imports prior sources read-only under
unique names. Native H components degree7/17/35 PASS1024shots each: independent
totalized model on split locators/all padding degrees, arbitrary paired
inverse, phase/scratch. Degree35 H forward41220T, paired82440T, Q897 including
852 publicbits, versusold50040/100080T. Reports
`results/fast_sylvester_H_checks.json`, `fast_sylvester_H_profiles.json`.
Degree35 fixtures all512 happened to compute trueH; degree7 has1 unsupported
interior gap, so this is not a repaired H acceptance claim.

Historical implementation design notes (NOW IMPLEMENTED and validated as
`repaired_sylvester.py`; latest state at top). Existing read-only source:
`../../5/dqi-core-phase/sylvester_circuit.py`. Need import with unique module
name and local outputs, no writes there. Exact source labels/operand order:
m.prod_rt(OUT,A,B). Do not misread triangular() target/source argument order.
Suggested native implementation: compute repair[j] = rlead_zero AND
rnext_nonzero AND no earlier repair bits; leave tags through phase. Conditional
rotate r right1. For c: rotate whole bank left1, then if k>0 swap c[k-1] with
c[-1]; gamma=c[k] is now at the last word. Cswap that with initially zero
shared gamma bank. Recurse. Restore r before correction; add
repair*gamma*x^k*r to low result. Undo source division as before.
Reverse H applies correction, recursive inverse, gamma/c swap, undo c rotate,
restore r, and uncomputes source-only repair flag. Need explicit inverse and
arbitrary-input paired tests, as well as natural/padded H correctness.

Avoid naive rotation of the whole2n-1 output bank for the x^k correction.
More efficient exact correction: rotate source r RIGHT by k (completion
count bits[1:]), multiply each rotated field coefficient by repair*gamma in
a12-bit temporary, XOR it into low_output[i]. For i<k the cyclic rotation
wrapped high coefficients: cancel at output[i] and add at output[i+n].
Both toggles cost one Toffoli/bit: CNOT(output[i],output[i+n]),
CCX([k>i],productbit,output[i]), CNOT again. [k>i] is completion
flag[35-2*i-2] at full degree35, or zero for a negative index. Undo r rotation.
This costs54nT for products+wrap, but source rotations cost about half the
whole-output rotations. Only emit k bits possible on repair: k<=n-2.
Mask gamma by current repair (12ANDs) so retained earlier gamma does not
affect later nodes. Source/payload conditional rotations and index selection
must be charged. k selection is linear in unary degree flags:
[k=h]=flag[35-2h] XOR flag[35-2h-2], with flag[35]=1 and negative flags0.
For payload suffix shift, conditional swap c[k-1] with last uses these bits
AND repair. Reverse loop order on inverse for arbitrary non-unary flags.

The H component/model and both full native jobs completed and the improvement
was promoted to best as recorded at the TOP. No GPU job remains from them.

`search_safe_schedule.py` all99 schedules saves only340T atQ1004; no rebuild.

Best is the accepted repaired-H package10.128965; first-crossed milestone10
is preserved. Milestones11–20 and score20 are not achieved. Continue.

## Latest post-acceptance research (2026-10-25)

`tensor_flip.cpp` exact field-multiplication decomposition search completed
59,395,257 flip/plus transformations in120seconds, seed2026102501. Best rank42
unchanged. `results/tensor_flip_search.json` includes source hash and scope.
No field arithmetic replacement is selected.

`chem_frontier_fit.py` fresh random R4/B240/C88, seed2026102503,12000 float32
steps, joint signed tensor/native Lambda/SOS-gap objective, COMPLETED170sec.
Final float63.5026mHa, rounded173.9854mHa, exact-word PSD-refit132.3009mHa,
refit lambda_eff78.6989. Fails physics and normalization; no native circuit.
`results/fresh_R4B240C88.json/.log`; checkpoint
`parameters/fresh_R4B240C88_last.npz` SHA256
a92cfede28c28d0d5a21a7daec661e9a9eb385c25c4d4bce8d7ecdda4573072c.
No unchanged restart is justified by this finite screen. GPU is idle.

NEW projective modular-inverse components COMPLETED and passed native tests.
`projective_polynomial_inverse.py`, `check_projective_polynomial_inverse.py`.
Latest (v2) degree3/7/17/35 all pass512 shots each for promised forward and
arbitrary paired inverse with non-unary flags, phase/scratch. Degree35 Q939,
T71782 forward/143564 paired, artifact
`candidates/projective2_poly_inverse_n35.kmb` SHA256
2939a88f239d6c3aff1c5aa3f2a9f35a7dd5ba74c9834a4f5c116a35008ba3aa.
Reports `results/projective2_poly_inverse_checks.json`, log
`results/check_projective2_poly_inverse.log`. Prior v1 Q939/T72202 was also
validated; source archived `results/code/projective_polynomial_inverse_v1_validated.py`
and report `results/projective_poly_inverse_v1_checks.json`.

Identity: F=c K(q),G=c K(tail q). Reverse q AND invert c:
Fnew=F/c^2 (same roots), Gnew=G^-1 mod F. No alternating head scales.
Retained inverses only at stable-source quotient division. Padded real-prefix
reversal uses partial right shifts then whole reversal. Both head suffixes
are zero except sole fake scalar in q1 bank; after shifts fake scalar is at
j=n-1-(n mod highest_power_of_two(m)). For n35 only3 possible positions.
Selectors are XOR-linear in unary degree flags. Park fake scalar in terminalc,
invertc, route it to canonical head[d].q1, inverse extract. Literal inverse
reverses all steps and handles arbitrary flags. Classical3239 n1..79/all d
permutation proof screen PASS: `results/projective_fake_scalar_shift.json`.
Interior Euclidean gaps remain unsupported; no complete submission claim.
T71782 is still120 above unchanged46-to11 scan break-even71662.5/direction.
Partial zero-prefix unpadding could save2424T whole lifetime, but any gap
repair must also be charged; no full build justified yet.

NEW `reciprocal_inplace.py` implements66T in-place GF4096 inversion instead
of102T low-space context+swap, by TRANSPORTING retained norm chain correctly:
N6 <-> Y6 (N^-1); M3 -> S3^2=M^-1; S3 -> M3^4=M^-3. All these are Clifford.
After source/reciprocal exchange, the transported chain is exactly the chain
of the new source, enabling measurement erasure with zero added Toffolis.
`check_reciprocal_inplace.py` exhausts all4096values/zero/two seeds, checks
involution and clean scratch/phase. Unconditional Q45/T66/8192shots PASS;
controlled Q46/T91/16384shots PASS. Reports
`results/reciprocal_inplace_checks.json`, log
`results/check_reciprocal_inplace.log`. Artifacts
`candidates/reciprocal_inplace{,_controlled}.kmb`.
This explicitly fixes chain state; do NOT simply enable unsafe retention at
source-swapping call sites. Old accepted source remains unchanged.

## Latest chemistry representation screens

`chem_sparse_pair_screen.py` completed40 fixed support screens, K52/104/152/
228/304 off-diagonal pairs per group,8 seeds each, plus shared76 coordinate
projectors. These dictionaries have exact3-bit Givens encodings, but disjoint
pair groups force cross-group target entries to zero. Best rigorous residual
floors are317.76/291.19/282.49/266.36/247.30mHa respectively, so all fail20.
`results/sparse_pair_capacity_screen.json/.log` (log is
`results/chem_sparse_pair_screen.log`). This excludes only saved disjoint
supports, not overlaps or denser low-bit dictionaries. Full target norm12.1842,
offdiag/offdiag targetnorm.84566; four-distinct-orbital coupling norm.34317.

COMPLETED numerical GPU session89026 under `/tmp/r6-gpu.lock`:
`chem_density_bootstrap.py --B240 --C88 --cycles10 --tag density_R4B240C88`,
log `results/density_R4B240C88.log`, final report sameprefix.json, checkpoint
`parameters/density_R4B240C88.npz`. New density/Schur factorization of fresh
signed native G: D=G[diag,diag]>0, first76 generators Xd=[sqrt(D),D^-1/2 C]
reproduce density and density/offdiag coupling exactly. Remaining276
rank-one modes of signed Schur complement fill4*C352. First block gets76
physical density generators +12 residual ones; remaining blocks88 each.
Build candidate pools from eigenvectors of16 high-strength matrices perblock,
coordinate vectors and immutable accepted719 atoms. Encode/decode8bits BEFORE
jointOMP selection (ridge.001 avoids near-dependence), then full PSD/rank88
refit10cycles and native gap projection. This is a different bootstrap from
random-gauge or commutator clustering. It completed12.5seconds: Schur target
truncation7.45237mHa, actual dictionary/refit591.083mHa, lambda_eff140.353.
Fails accuracy and cost; no circuit or native acceptance. No unchanged restart.

`chem_condition_grow.py` COMPLETED52.3seconds, no GPU job remains. Grows best
original-tree R2/B479/C175 native-physical checkpoint toB511 with32 exact atoms
perblock in8batches, pricing ridge residuals WITHOUT outside-span projection
so cancellation-repair directions are eligible. Duplicate projectors excluded.
Seven declared ridge levels0..1e-3, each fixed-generator fit and4+16 PSDcycles.
Best error22.1378357mHa (lambda_eff96.4491), lowestlambda71.4532 at26.7321mHa.
Initial22.5122444mHa/lambda71.4770. No passing or cost-improving model.
`results/condition_R2B511.json/.log`, checkpoints `parameters/condition_R2B511*`.
Fixed-generator exact refits cannot redistribute weights over linearly
independent appended atoms: their exact new coefficients are zero. Finite
ridge relaxation changes factors and can affect weights but did not solve
this screen. Do not confuse L2-conditioning gains with physical L1 gains.

## Open algebra idea, not implemented or profiled

Jacobi-to-unnormalized-continuant synthesis may avoid long scalar scales.
For f=a0/(1+b0*x+x^2*a1/(1+b1*x+...)), set t=1/x and R=f(1/t)/t.
Choose k0=1/a0 and ki=1/(ai*k_(i-1)), qi=ki*(t+bi). Then
R=1/(q0+1/(q1+...)). Transform each (ai,bi) into(ki,ki*bi), arrange
as degree-one quotient heads with terminalc=1, then inverse the unnormalized
polynomial EEA extractor to produce (F,G), with R=G/F. Reverse coefficient
lists to get the ordinary locator/evaluator, allowing their common scale.
One extra12-bit explicit leading/terminal word suffices dimensionally.
Need complete coefficient conversion AND literal inverse costs; a source
reciprocal swap may share the new chain-transport primitive. This is not
the prior symbolic short-block `continuant_coordinate_profiles` screen.
No width/T or score claim; zero charts still require explicit treatment.
The phase-only original mask interface and gap statistics remain mandatory.

No processes are live from this post-H turn. Best remains10.128965274112238;
first-crossed10 archive preserved; goal20 remains ACTIVE,11–20 missing.

## Latest wave: completed reciprocal and quotient components, positive chemistry

Accepted best and immutable score-10 milestone are UNCHANGED. No new full
submission verification was warranted. Goal20 remains ACTIVE;11–20 missing.

`copied_prefix_inverse.py` implements F^-1=F*(F mod x^ceil(n/2))^-2
mod x^n by copying/inverting a short prefix, applying its squared unit-series
multiplier, and erasing the copy against the output prefix. Ignored physical
F0 remains arbitrary. `()` must enter known_zero BEFORE recursive inversion.
`check_copied_prefix_inverse.py` native reciprocal/phase/scratch/involution
PASS at n6/14/34/70, generic/zero F1, one-copy/recursive,512 shots each.
N70 one-copy genericQ1251/T22764, zeroF1Q1239/T20286; recursive generic
Q1623/T20412, zeroF1Q1563/T16590. Report
`results/copied_prefix_inverse_checks.json`. This trades width for arithmetic.

IMPORTANT profiling correction: `results/copied_core_profiles.json` only
changed the outer arithmetic object. Inner `fused_gap_block` and `gap_serial`
constructors still used the prior class. Its T790980 is the earlier616T parity
saving, NOT a reciprocal-copy gain. Archived source
`results/code/profile_copied_core_outer_only.py`. Corrected
`profile_copied_core.py` overrides the nested constructors too. Report
`results/copied_core_profiles_nested.json`: limits940/980/1004/1080/1200 give
(Q,T)=(1004,796740)/(1004,796308)/(1011,794364)/(1017,794052)/(1017,794052),
versus selected core791596. New copy used48 or88 times at higher budgets,
but this small-head realization LOSES. Do not build a full candidate from it.

`quotient_jacobi_synthesis.py` implements the prior algebra sketch. k0=1/a0,
ki=1/(ai*kprev), qi=ki*(t+bi), reconstruct R=G/F by inverse unnormalized
polynomial EEA. Generic n35Q898/T42836 one direction, versus prior monic
Jacobi synthesis53822. All n3/7/17/35 native generic rational recurrence,
paired arbitrary-zero restoration, phase and scratch PASS512 shots each.
`results/quotient_jacobi_synthesis_checks.json`. Optional `before_swap`
callback added to `reciprocal_inplace.apply` shares chain inversion with
bi scaling. Default None path emits the previous operations unchanged.
The earlier post_H manifest still binds the EARLIER source; do not pretend
its reciprocal source hash matches the new optional-callback source.

`completed_quotient_jacobi.py` now handles all COMPLETION degrees, while
interior Euclidean gaps remain unsupported. Store35 zero-a flags, c=1 for
full degree or0 for completion, and arm ONLY the first zero quotient as t.
Totalized EEA reconstruction produces F=t^(n-d)*Ftrue,G=t^(n-d)*Gtrue.
The n35 native component is Q933/T43255 (paired86510),512 shots/circuit,
SHA475d13635ed390a4395f2e65d1e773373731cb883808f959742c8d557df88a10.
All n3/7/17/35 PASS: independent Jacobi recurrence cross multiplication,
all completion degrees, exact common monomial padding, arbitrary heads on
paired inverse, phase and scratch. `results/completed_quotient_jacobi_checks.json`.
This is a component, NOT a complete DQI oracle. With direct Forney G/F',
common t^m cancels at nonzero roots; inverse-cofactor/H endpoints still need
correct unpadding. Degree-two quotient/3-parameter gap packets and deferred
second-pivot gaps remain potential future work, not implemented here.

### New PSD-generator chemistry representation

`chem_positive_generator_screen.py` necessary positive-map test found no
obstruction:332 unit tests, min output eig+.00302192, no violated pair
inequality. `results/positive_generator_screen.json`. This test alone is
not feasibility, but subsequent ADMM DID construct a positive factor.

`chem_psd_generator_rotation.py`: fixed-tensor polar rotations plus spectral
PSD targets. Accepted source C350/700 after200 iterations still has negative
Frobenius .030025/.028147. Fresh native best rank352 PSD source error is
7.3611395565mHa;500 alternating projection steps leave negative norm.014518,
and its PSD projection error40.0097mHa. Finite rotation screens, not cone
infeasibility proofs. Reports `results/psd_generator_rotation.json` (original
incumbent run), `results/psd_generator_rotation_native.json`.

`chem_psd_generator_admm.py` accelerates the fresh native source using dual
updates: PSD(Y+dual), polar fixed-Gram projection of Z-dual. At step190 all
352 generator matrices are PSD. Step200 minimum eig+1.55536212e-7, negative
norm0, native tensor error7.3611395565mHa, fixed-Gram drift checked<1e-8.
Report `results/psd_generator_admm_native.json`, parameter
`parameters/psd_generator_admm_native_C352.npz`. This establishes a numerical
PSD-generator factor under the fixed native shift, independent of a finite
rank-one dictionary. It is not a submission or a resource improvement.

`chem_positive_factor_physics.py` intrinsic trace normalization: for any
nonnegative projector coefficients, column L1 equals Tr(Xc). Convex gap
LBFGS400 at this positive factor gives minimum gap4.08171097, Lambda283.22309,
lambda_eff47.9104291; two-body7.361mHa. Starting negative trace offset had
gap174.65, so offsets must be optimized. Report
`results/positive_factor_physics.json`. No finite dictionary or oracle yet.

`chem_positive_dictionary_fit.py` exact native tensor and eight-bit metrics:
R4/B240/C88, common-mean eigenvectors plus generator spectral atoms,
nonnegative NNLS initialization.4000 fixed-generator fit steps finish
521.372mHa float /517.604mHa exact words (fails).10000 tensor-loss steps
with log weights improve64.0097mHa /160.365mHa. Direct projected nonnegative
weights allow activation at zero;12000 steps finish38.0827mHa /152.099mHa.
This is a material optimizer change, not a rerun with a different random seed.
Reports `results/positive_{R4B240,tensor_R4B240,projected_R4B240}.json`.
`chem_positive_lbfgsb.py`500 bound-constrained native float64 iterations
further reduce float35.9499mHa; rounded171.092mHa. Still fails. Report
`results/positive_lbfgsb_R4B240.json`.

`chem_positive_dictionary_refit.py` checks best projected integer geometry:
union-span necessary floor92.3555mHa already rules out20mHa for THESE words.
Signed rank88/block PSD fits1+4+20 cycles reach121.863mHa; fullC240/block
relaxation reaches118.684mHa. `results/positive_projected_refit.json`.
These are fixed-word bounds and finite fits, not shape-wide impossibility.
The positive source exists, but this compact quantized dictionary fails.

`chem_coupled_word_repair.py` screens all75-choose-2 angle pairs and4 signs
using full coupled Gauss-Newton matrices, then exact row residual changes
and full native physics. Exact single-row prediction agrees<1e-9. On B703
low-normalization checkpoint, cap51 yields no feasible move; slightly relaxed
cap51.05 allows271 rows and reduces20.144103→20.141605mHa before hitting
normalization boundary again. No pass or resource improvement. Reports
`results/coupled_B703{,_relaxed}.json`. No unchanged restart justified.

CURRENT RUN: native integer refinement of the NEW positive-source dictionary,
not an old repaired geometry. `chem_positive_integer_refine.py` is a local
copy of read-only round5 chem-audit refiner with HERE redirected locally;
all outputs stay agent20. Exact native signed PSD rank88 block refits and
1/2/4/8-angle proposals, finite300 steps, unbounded physics for capacity
research. Init `parameters/positive_projected_refit_ranked.npz`, output
`parameters/positive_native_integer_R4B240.npz`, log
`results/positive_native_integer_R4B240.log`, session93006, GPUlock.
Inspect process/log/report before launching anything else. Native full metric
has NOT run and no candidate was promoted. All other sessions in this wave
completed successfully; GPU is used only by this finite refiner.

### Latest generalized-quotient classical model (new, quantum circuit unbuilt)

`generalized_quotient_jacobi_model.cpp` eliminates the fixed first/second
zero-pivot distinction at the mathematical level. For valuation v and L=v+1,
write f=x^v*a/(1+b1*x+...+bL*x^L+x^(L+1)*g). The packet consumes2L moments
and has v imposed zeros plus(a,b1,...,bL), exactly2L physical head slots.
For R=f(1/t)/t this is R=a/(t^L+b1*t^(L-1)+...+bL+Rnext).
Convert with ki=1/(ai*kprev), independent of quotient degree, then use
ordinary polynomial continuants. Reversing ALL variable-degree quotients
keeps F exactly and gives G*Greverse=1 mod F, including interior gaps.
For the generated syndrome convention, e=G(a)/(a*F'(a)), and therefore also
e=1/(a*Greverse(a)*F'(a)). This was checked directly at every actual root.

Original source archived at
`results/code/generalized_quotient_jacobi_model_law.cpp`:100000 predeclared-law
vectors plus256 forced leading valuations1..4, ZERO wrong vectors,991 vectors
with gap, maximum quotient degree5. Report
`results/generalized_quotient_jacobi_model.json`.
Current source adds explicit completion-degree fixtures d0..35; separately
run with0 law shots,288 degree fixtures plus256 forced valuations, ZERO wrong.
Report `results/generalized_quotient_jacobi_degree_fixtures.json`.
Both check an independent generalized Jacobi polynomial recurrence, exact
F symmetry under quotient reversal, the original root/magnitude vectors,
and the reversed-cofactor magnitude formula. These are classical model
checks only: no quantum resources or submission acceptance are claimed.

Potential next implementation, not yet coded: a shared degree1/degree2
polynomial quotient peel. Let input F degree n,G nominal n-1 and flag=G[-1]==0.
The common low remainder is F[:-3], length n-2. Its entire tail convolution
can use the SAME G[:n-2] on both charts and a selected three-word quotient:
normal(q0,q1,0), gap(q0,q1,q2). Only the normal boundary F[-3] needs a separate
remainder update; on the gap chart it stores q0. Use pivot
L=G[-1]+flag*(G[-1]+G[-2]), and corresponding selected sublead for the top
quotient recurrence. Masking these few coefficients avoids duplicating the
long polynomial product. Lower-degree physical layouts, cross-block deferral
and completion routing must still be implemented and charged. A degree2
packet at the last slot of a fixed block must wait/carry into the next block;
do not silently balance it as2/2. This is a possible bridge from the now-valid
all-gap model to the cheaper quotient-synthesis components.

Latest running chemistry status: session93006 at step200/300 has exact-word
error58.2268mHa,4242 accepted row changes. Native gap-projected checkpoint
metrics are Lambda394.475,gap4.12146,lambda_eff56.8739. It still fails20mHa.
The unprojected main-offset metrics in intermediate log rows are deliberately
poor because this run explores tensor capacity; the `metrics.refit` fields
include native gap projection. Let the predeclared finite run finish, inspect
its final report, and do not restart it unchanged. Parameter checkpoint is
mutable until terminal state. No other jobs are live from this wave.

Completed immutable research artifacts and source/byte hashes are recorded in
`results/reciprocal_quotient_positive_research_manifest.json`. Accepted best
and score10 milestone hashes were rechecked unchanged. Goal20 remains active.

## Latest head/quotient wave (supersedes the running status above)

Accepted best STILL10.128965274112238, independent10.144446072819372.
No new full-host candidate accepted; milestones11–20 remain missing.
The prior chemistry integer run finished300 steps at56.7293350141mHa;
gap-refit56.7284124911mHa/lambda_eff56.860314. Fails20mHa; no unchanged restart.

Completed quantum components, all with native phase/scratch/literal-inverse
checks and independent scalar arithmetic. These are NOT submission scores:

* `shared_gap_quotient_peel.py`: one shared long remainder product for degree1/2.
  n3/7/17/35 pass, n35 Q947/T1972 per direction. Original source archived
  `results/code/shared_gap_quotient_peel_v1.py`. Current source routes the third
  coefficient through a zero word; it also passes n3/7/35 at the SAME peak/T.
  New KMBs are `candidates/routed_shared_gap_quotient_n*.kmb`. Original paths
  were rebuilt from v1 and SHA-checked against the original report.
* `continuant_frame_affine.py`/`continuant_head_frame.py`: compact A/D/p and
  direct affine formulas; native parameter-XOR tests pass. Variants archived
  under results/code. Current b2 single N70 block Q935/T43644; larger symbolic
  formulas lose. SMALL N14/M127 oracle Q245/T37660 passes4096 shots on its
  model/restoration contract,153 truth errors on unsupported gaps. Cache b4
  capacity4 sweep was explicitly aborted, not a resource lower bound.
* `compact_lookahead_parameters.py`: leading-zero two-pair parameters in head.
  N6/14/34 components pass. SMALL N14 Q243/T51636 passes4096 shots, zero truth
  errors including283 forced unsupported labels. Full N70 extraction401370T
  per direction; old unshared synthesis207629T. ProfileQ915 excludes the scan.
* `counter_quotient_jacobi.py`:35 completion flags become6 count+1 armed bits.
  n3/7/17/35 pass all completions, independent recurrence and arbitrary paired
  nonunary zeros. n35 Q905/T44206. Its convenience model is only unary-valid;
  `prefix_counter_reference.synthesize` implements the totalized scalar map.
* `high_copy_jacobi.py`: retain upper b raw head moments, compute affine P
  as high coefficients of A[1:]*raw_high, share long tail; restore head briefly
  to erase copy. Native b3 N8/14/34 pass. Single N70 b2/3/4/5/6 qubits
  935/960/967/992/1017, Toffolis43590/51952/58754/66030/73726.

Narrow prefix line: `prefix_high_copy_jacobi.py` uses a simple early cofactor
frame and changes the raw-high bank to P in place. Zero A_b is armed, so it
always restores but may decode incorrectly. Native b3 N8/14/34 pass. Later
blocks retain full completion charts. `profile_prefix_counter.py` b6/t18/s2
Q958/projected completeT911108 (80*4095 scan); inverse limit935, t20/s3 gives
Q942/T915560. Corrected100000-vector census:785 wrong for b6/t18/s2,813 for
b6/t20/s3. `prefix_counter_model_b6_t18.json` is OBSOLETE (wrong late frame);
the corrected report is `_v2.json`. No full native benchmark justified.

`prefix_counter_oracle.py` is a complete research builder, scanning roots
a=alpha^(+j), e=G(a)/(a F'(a)), virtual F1/F0/G0 accumulators, dual classical
1/a masks. Terminal/count words remain reserved. Initial twice-phase failure
was private classical measurement garbage; current builder clears private
scan bits. SMALL N14 b3/t3 Q263/T39476 passes2048 shots: all input/mask/model,
scratch/private bits and twice-phase checks zero.271 unsupported truth errors.
SHA b488a100e2515be28f34649e20f7487dfb14676f44cc8a0cdbe586f5bd785947.

NEW `mixed_quotient_jacobi.py`: normal two-linear versus one-quadratic packet
reconstruction, with shared long products. Convert to quotients, retain one
gap flag per suffix packet, count completion excluding real gap markers.
After the shared normal inner fold, swap tail banks on gap, then the shared
outer fold. Final degree2 uses two known-zero low words. Native full-mixed
n3/7, even n2/4/6, n17/35 with18 normal prefix slots all pass variable-degree
rational recurrence, all completion degrees/common padding, arbitrary paired
zero states, phase and scratch. Reports `mixed_quotient_jacobi_{odd,even}_checks`
and `mixed_quotient_jacobi_p18_checks`. n35/p18 Q961/T56588 per direction with
inverse policy1004. Selected reciprocal must use `reciprocal_inplace.apply`'s
explicit control while norm setup/erasure stays unconditional. Passing a
selected emitter to the whole reciprocal hits unsupported conditional CCZ.

`mixed_prefix_oracle.py`: narrow b6 prefix through18 plus compact leading-gap
two-slot suffix and mixed synthesis. SMALL N14 b2/t2 Q280/T47872 passes2048
shots incl293 arbitrary syndromes;168 truth errors, none supported. Routed
variant Q280/T46432 passes the SAME cases, SHA
0f086c09c8d78d609d3cf19e1367497e50c767f3a2ee4e18f57b0123784ea7e6.
Full arithmetic profile with inverse limit935: Q943, core639936T, projected
complete967536T. C++ `mixed_prefix_model.cpp` finds658 wrong vectors/100000,
seed2026102621; this still prices above the accepted DQI. No full native build.
`check_mixed_prefix_model.py` compares every coefficient with the scalar model
on every completion degree N12/N14/N70, arbitrary and nonunary-zero inputs.

`routed_quantum_control.py`: selected field scale through one routed zero word
instead of two copied factor words. Target/factor aliases fall back to old
implementation. Install in a NEW process only; patches imported references,
never prior files. Mixed saves5760T, unchanged peakQ943. Accepted H saves only
48T, unchangedQ1004. `profile_routed_H.py` wrote identical generated H source;
all accepted source hashes were rechecked. No rebuild for that negligible gain.
Broader exact quadratic searches completed: mixed-operand field multiply
stays42 after31,272,667 transformations; joint Forney norm/pairing stays45
after22,359,078. Finite searches, not lower bounds. No arithmetic replacement.

Chemistry ADMM: v1 fixed rank352 tensor-factor Procrustes plus exact-word
least-squares projection/local row moves finished100 steps at387.9mHa,
minimum reported93.0mHa, worse than56.73mHa init. Weighted v2 uses right metric
G^(1/2)+.005I and dual step.1:150 steps finish267.6mHa, best stays56.729335mHa.
Sources archived `results/code/chem_native_factor_admm_v{1,2_weighted}.py`.

The following global-atom run is now COMPLETED (session82855), finite150-step
`chem_native_factor_admm.py --dual-step .1 --global-atoms --tag
native_global_atom_admm_R4B240`, GPUlock, log of same stem under results.
It finished at409.414904mHa; best stays56.729335mHa init. No unchanged restart.
New global K-SVD atom
proposals use the lowest eigenvector of unpack(w_b*residual-||w_b||²P_b), the
exact continuous atom minimizer at fixed coefficients for unweighted factor
loss. Round whole-vector mixtures at.25/.5/1 and accept only exact-word full
refit improvements. This is distinct from another unchanged local word run.

Possible next quantum step, UNBUILT: triple packets support gaps at slot0 or1,
leaving slot2 for the next block. Heads: normal(a,b,c,d,e,f), first-gap
(0,u,h,e,a,b), second-gap(a,b,0,u,h,e). First-gap cofactor is
A=t³+(e+b)t²+(h+be)t+(bh+a), D=t²+et+h, p=u*armed(a).
Second-gap has constant bh+u instead, D=t+b padded with zero second
x-coefficient, p=a*armed(u), and tail reciprocal shift3 instead of2. These
fit six head words: first gap compute p in u, A3 in a, A2 in original zero,
A1 in b, then route slots(0,1,4) to normal layout; second gap compute p in a,
A3 in u, A2 in h, A1 in e, then swap slots1/3. Compute p before changing its
source. First/second head charts choose leading-gap+single versus
single+leading-gap on six words. Native completion/inverse tests still needed.

UPDATE: `triple_gap_frame.py` and `check_triple_gap_frame.py` now implement
the parameter frame just described. Native arbitrary-target parameter XOR
passes1024 regular/gap/completion shots plus512 arbitrary paired marker and
flag states, all source/output/phase/scratch checks. Q244/T4246 includes72
input and72 output qubits plus2 public chart bits. Circuit SHA
89229483b5497b5b536f18b3510afd3ff91db395b3fea17efad70a3d1231497f.
Head extraction, long-tail integration and three-slot quotient synthesis
remain UNBUILT. Only the cofactor frame is validated. Final source and byte
hashes for this wave will be in `results/head_quotient_research_manifest.json`.
No GPU/full verifier job remains after the terminal checks above.

## Completed generalized quotient wave (supersedes UNBUILT notes above)

Accepted best and milestone10 are unchanged: 10.128965274112238 default,
10.144446072819372 independent combination. No new full-host verifier run or
promotion in this wave. Milestones11–20 remain missing; goal20 stays active.
All native component and GPU jobs below have finished.

`triple_gap_head.py` now implements six-word head extraction. Selected v1
Q181/T34344 passed but lost badly; source/native bytes are archived with
`selected_v1` names. Shared encoded v2 Q145/T2432 is also archived.
Current canonical v3 converts the normal h/e tail back to ordinary Jacobi
coordinates, Q145/T2498 each direction, native1024 shots/direction all phase,
scratch and independent model checks pass. Forward SHA99f6a126589cef73dc96e2ba5a3ab7b0606c65f6950295904b4beab6d15d52af;
inverse SHA78338c9cd1a77a0e003aa0034a2b1a6cb330280496ee98bf195cceeac5dfeab7.

`triple_gap_block.py` combines this head with the validated triple cofactor
frame and variable tail shift. A dirty chart flag originally caused phase
errors: cleanup must erase the second flag before the first. Only the final
`check_triple_gap_block_fixed_v2.log`/checks report describes the corrected
PASS. N8/14/34: Q214/286/526, T13368/16388/32130 per direction,1024 native
shots each direction. N70 single block Q958/T66454 (inverse limit935), or
Q976/T65482 (limit1004). Full extraction Q958/T425579 per direction, before
synthesis/scan/inverse. Too expensive; gaps at each third block slot still
unsupported. The old phase probe is diagnostic evidence of the fixed bug.

`generalized_gap_quotient.py` reconstructs variable-degree1/2 continuants in
fixed slots, including arbitrary adjacent flags with skip priority. Native
v1 n35 Q989/T52454 per direction passed256 shots; original source archived
under results/code. Direct n2 terminal v2 avoids24Q. Current v3 shortens
sublead lifetimes and uses one aggregate product for the third shear:
n2/3/7 Q134/159/259 and T554/1166/4466, all256 native shots/direction plus
arbitrary paired flags/zeros pass. Sparse flags n7 T4340 also pass.
Large v3 is profiled, NOT native-tested: full-triple footprint Q954/T48978
for reconstruction only. Six-head footprint with all840 public bits and72
raw-copy bits Q981/T3636 at inverse limit1004. Not a full decoder.

`quotient_cofactor_frame.py` reverses generalized quotient packets and uses
normalized reversed continuants to obtain compact A/D/p. It retains one
gauge word, avoiding a resultant computation. Native arbitrary-target XOR
n2/4/6: Q181/278/375, T4120/9876/16100,512 regular/gap plus256 arbitrary
paired shots each, all phase/scratch/source/output checks pass. Arbitrary
zeros give a total inverse; completion decoding correctness is not claimed.
Report `results/quotient_cofactor_checks.json`.

`prefix_lookahead_quotient_block.py` integrates that frame with b2 lookahead
heads and an upper-triangular raw-head boundary map. Native b6 N14/34 passes
1024 shots/direction against an independent totalized scalar map:
Q328/T45258 and Q568/T62664. Reported `serial_matches:2` compares different
head encodings; it does NOT isolate tail correctness. A separate supported
serial-tail comparison remains useful. Single N70 b6 at inverse limit1024
Q1000/T102376, worse than old generic prefix Q958/T72470; no full build.

Chemistry controls are terminal. Current `chem_native_factor_admm.py` is v4;
global150-step v3 is archived as `results/code/chem_native_factor_admm_v3_global.py`.
The historical manifest's mapping of that report to the current source is
stale and must be corrected before rerunning its recorder. Weighted global
ADMM v4 with20 inner projection steps finished400 steps at201.517906mHa,
best still56.7293350141mHa initial. `native_global_atom_weighted_inner20`
JSON/log/checkpoints retain this failed control. No unchanged restart.

`chem_positive_dictionary_fit.py` now permits B<76 initialization; original
source archived `results/code/chem_positive_dictionary_fit_v1.py`.
R11/B64/C32 finished6000 steps: float174.296616mHa, exact-word NNLS267.217732mHa.
`chem_fixed_dictionary_union_floor.py` gives the saved exact-word dictionary
an unrestricted signed-weight union-span floor158.634857mHa. Refitting those
vectors cannot pass20mHa; this is not a shape-wide exclusion.
The fixed PSD factor's trace/Frobenius ratio is9.80903129, a necessary common
generator-rank bound, not a feasibility guarantee for any R/B/C shape.

Current-wave bindings are in `results/generalized_quotient_research_manifest.json`.
The older head-wave manifest describes a historical snapshot, not all later
source revisions. No accepted bytes were changed by this research.
