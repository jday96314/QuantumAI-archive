# Adaptive dyadic inverse square root — round 5

This experiment has not found the requested substantial cost reduction. The
best shared-coordinate construction reproduces the incumbent's 90-qubit,
1793-static-Toffoli resource point. Variable-depth local coordinates are
functional but cost more. The 165,000, 155,000, and 140,000 targets are not met.
The verified round-4 incumbent is retained unchanged, including a local reference
copy with SHA256
`27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1`.

The reusable improvement is [stream_erase.py](stream_erase.py): streaming one
ANF bank reduces lookup-erasure workspace without increasing its Toffoli count.
For an 11-bit address it needs 31 rather than 83 scratch wires and still costs
83 Toffolis. This helps the losing linear and mixed circuits; it does not lower
the incumbent's peak.

## Final native verification

<!-- NATIVE_RESULTS -->
Both frozen experimental candidates **PASS all 200,540,928 native shots each**.

| Circuit | Q | Mean Toffolis | Rotations | Full task cost | Worst error / tolerance |
| --- | ---: | ---: | ---: | ---: | ---: |
| Retained round-4 incumbent | 90 | 1782.000090874 | 34 | 182956.52716064 | 0.979671439549 |
| Shared prefix refinement | 90 | 1782.000221551 | 34 | 182956.53892153 | 0.979671439549 |
| Variable-depth local coordinate | 92 | 1805.999167422 | 34 | 188762.16929350 | 0.938422588735 |

Native evidence: [final_shared_verify.json](results/final_shared_verify.json),
[final_local_verify.json](results/final_local_verify.json), and the copied
[incumbent verifier result](reference/incumbent_verify.json).
The [audit](results/audit.json) binds both complete results to candidate bytes,
checks exact shot counts and numerical synthesis margins, confirms identical
rebuilds, and confirms the incumbent is unchanged. Total final native shots:
**401,081,856**. The shared circuit is approximately 0.01176 task-cost units
above the incumbent because of finite measurement-branch averaging; this is
not an architecture improvement.

Fresh C++ refits reproduce both final leaf lists exactly; see
[fit_reproduction.json](results/fit_reproduction.json). Independently rebuilt
complete circuits are byte-identical to both frozen native-verified candidates.

Frozen candidates and SHA256:

- [final_shared.kmb](candidates/final_shared.kmb):
  `a8c1b24b1fad8c5a63a2f4c7f20e79b7692b051290d9cbe784ac08a03801be40`.
- [final_local.kmb](candidates/final_local.kmb):
  `2305f7c24debd5217efe5a6813bfe3e86cf63e5afca94c2af4b00f60f5c864f8`.
<!-- END_NATIVE_RESULTS -->

Only the two `final_*.kmb` artifacts are final experimental candidates. Other
compiled resource points below are exploratory: their numerical models check
all inputs, and selected complete circuits have shortened native checks. They
are not being represented as full native passes.

## Scope and independence

All writes are confined to `initial_experiments/5/inv-dyadic`. No host verifier,
previous round, completed round-5 experiment, or submission archive was modified.
The selected source chain and numerical/workspace experiments in round 4 were
read before implementation. No other round-5 experiment's contents were used.

[fit.cpp](fit.cpp) generates coefficients without reading the incumbent's table.
[certify.py](certify.py) independently computes integer Horner arithmetic and
the exact fixed-point target. [build.py](build.py) is a new builder; it does not
call or rewrite the incumbent's whole-circuit builder. It reuses the established
normalization, radix-four multiply-add, unequal-width adder, and ANF-erasure
components through read-only imports. Final source and dependency hashes are in
[source_manifest.json](results/source_manifest.json).

One relevant baseline finding: the incumbent's nominally 512-row coefficient
table already compresses to **269 binary-prefix leaves**. Its lookup is not a
512-Toffoli flat lookup. The independently fitted 269-leaf shared-coordinate
table differs from the incumbent in only one expanded row (constant and slope).
The comparison was performed after fitting; the incumbent table was not a seed.
See [incumbent_table_comparison.json](results/incumbent_table_comparison.json).

## Integer approximation and dyadic partitions

For input `x`, define

```
k = 21 - floor(log2(x))
M = x << k                       # 2^21 <= M < 2^22
m = M / 2^21
p = k & 1
s = k // 2
g_p(m) = (1 if p else 1/sqrt(2)) / sqrt(m)
```

A leaf record in `models/*.json` is

```
[parity, lo, depth, coordinate_origin_depth, degree, C0, C1, C2, fit_error]
```

It covers `M in [2^21+lo, 2^21+lo+2^(21-depth))`. Thus every boundary is exactly
a binary prefix. The complete leaf lists, not a floating boundary convention,
define the partitions. Local architecture leaves use their own midpoint;
shared architecture leaves use the midpoint of their coarse ancestor.

For coordinate depth `f`, the signed integer coordinate is

```
z = (M mod 2^(21-f)) - 2^(20-f)
```

It remains in input wires. The multiplier uses its two's-complement encoding,
including sign extension in its final radix-four digit. Define `H(A,z,e)` by
recoding `z` into digits `d_g in {-2,-1,0,1,2}` and adding, for every nonzero digit,

```
sign(d_g) * round_half_up(A * 2^(e + 2*g + log2(abs(d_g))))
```

Every row rounds individually. Rounding an aggregated real product is a
different approximation and was not used for scoring these models. The integer
quadratic is

```
B = C1 - H(C2,z,w1-w2-21)
P = C0 - H(B,z,w0-w1-21)
phase = (P << s) / 2^(w0+10)      # turns, modulo 1
```

Linear leaves set `C2=0`. Mixed tables combine both degrees in the same schedule.

The fitter seeds derivatives at leaf midpoints, shifts them to the selected
coordinate origin, and quantizes each coefficient. It searches integer slope
offsets of radius 18–48 and curvature offsets `-3..3`; it optimizes the constant
against the extrema of the actual rounded integer outputs. A 192-point screen
ranks candidates, the best eight are checked on every valid normalized integer,
and a failing leaf splits in half. The screen can miss a better fit; no claim of
global coefficient optimality is made. Derivative/curvature seeds and sampled
integer errors generate the partition, and exhaustive integer error refines it.

Parity-dependent reachable normalized integers are enforced: odd `k` requires
even `M`; even `k` permits `M <= max_x` or multiples of four. The subsequent
certificate checks every original `x=1..3,133,452` directly. Fit targets, depth
caps, search radii, and reproduction inputs are in
[specifications.json](specifications.json).

### Approximation frontier

There are 20 models. Selected points are below; the full machine-readable
frontier is [approximation_frontier.csv](results/approximation_frontier.csv).
`ratio` is exact phase error divided by the host's integer tolerance. PASS here
means the numerical model, not native circuit verification.

| Model | Fractional bits `(w0,w1,w2)` | Leaves / depths | Worst ratio | Model |
| --- | --- | --- | ---: | --- |
| Local quadratic | 29,19,9 | 163 / 6–7 | 0.938443 | PASS |
| Local quadratic | 28,19,9 | 184 / 6–7 | 0.938197 | PASS |
| More accurate local quadratic | 28,19,9 | 203 / 6–7 | 0.833801 | PASS |
| Local quadratic, narrower derivatives | 28,18,8 | 251 / 6–7 | 0.935044 | PASS |
| Local quadratic | 27,18,8 | 285 / 6–8 | 0.938423 | PASS |
| Local quadratic starting at depth 7 | 27,18,8 | 287 / 7–8 | 0.938423 | PASS |
| Shared coordinate at depth 6 | 27,18,8 | 518 / 6–9 | 1.340779 | FAIL |
| Shared depth 6, wider derivatives | 27,19,9 | 217 / 6–8 | 0.939299 | PASS |
| Shared depth 7, tighter margin | 27,18,8 | 287 / 7–8 | 0.938423 | PASS |
| Shared depth 7, fit target `5.9e-8` | 27,18,8 | 269 / 7–8 | 0.979671 | PASS |
| Shared depth 8, narrow derivatives | 27,17,7 | 512 / 8 | 0.930027 | PASS |
| Shared depth 7, narrow derivatives | 27,17,7 | 645 / 7–10 | 0.961598 | PASS* |
| Local linear | 29,20,10 | 1876 / 9–10 | 0.938167 | PASS |
| Mixed: final 1/32 of mantissa range linear | 28,19,9 | 213 / 6–10 | 0.9382 | PASS |
| Mixed: final 1/8 linear | 28,19,9 | 321 / 6–10 | 0.9382 | PASS |
| Mixed: final 1/2 linear | 28,19,9 | 1024 / 6–10 | 0.938206 | PASS |

`*` The 645-leaf model misses the tighter `5.6e-8` fitting target but passes
the actual `2^-24` numerical requirement. Three output-26-bit searches fail the
actual requirement, with ratios 1.112301, 1.062407, and 1.050795, despite up to
603 leaves. These are failures of the searched quantizations and rounding
schedule, not impossibility results for 26-bit approximations or dyadic circuits.

The 184-leaf local partition has 72 depth-6 and 112 depth-7 leaves. The final
shared partition has 243 depth-7 and 26 depth-8 leaves. The final variable-local
partition has 225 depth-7 and 62 depth-8 leaves. Thus asymmetric refinement is
actually used, rather than merely relabeling a uniform partition.

### Exact accuracy, overflow, and synthesis margins

The reference integer is `r=isqrt(floor(2^130/x))`; the target is `(r+1)//2`
modulo `2^64`, and tolerance is `r>>25`. The computed phase is an exact dyadic
integer in the same units. The certificate subtracts modulo `2^64`, checks all
3,133,452 inputs, and takes the **minimum absolute leftover tolerance over the
entire domain**, separately from the worst relative-error ratio.

For each final candidate the host synthesis rule is used unchanged:

```
eps = pi * min_margin_units / 2^64 / rotation_count
T_synth = rotation_count * (0.56*log2(1/eps) + 5.3)
work = mean_Toffolis + (T_synth + T_base)/4
```

| Model | Minimum margin (turns) | Rotations | Per-rotation eps | Synthesis + T-base work |
| --- | ---: | ---: | ---: | ---: |
| Final shared | 1.0839411508124785e-12 | 34 | 1.0015592812694554e-13 | 250.850211 |
| Final local | 2.2734988804756773e-12 | 34 | 2.1007080531902962e-13 | 245.763542 |
| Local 184-leaf | 2.210051911442912e-12 | 35 | 1.9837379568688948e-13 | 253.389534 |
| More accurate local, 203 leaves | 6.28724546383741e-12 | 35 | 5.64341833157353e-13 | 245.998630 |

The additional accuracy of the 203-leaf model saves only about 7.39 units of
synthesis work relative to the 184-leaf model, while adding 19 leaves. Nominal
relative accuracy therefore does not predict total cost.

Every certificate also records intermediate slope and mantissa ranges, rounded
product magnitudes, product-sign checks, and signed-accumulator bounds. Final
local and shared models fit signed accumulator splits `(18,9)`. Their products
have magnitudes at most `(261376,380)`, below `(2^18,2^9)`, with no sign
violations. Mantissas stay in `[67108897,134217726]`, making the aliased top
input bit exactly one before phase shifting. Native simulation checks the
modular implementation and scratch rather than relying on these ranges alone.

An actual overflow hazard was found: for depth-6 leaves, a nonnegative loaded
constant does **not** imply that the affine residual remains nonnegative after
Horner. The 184-leaf model reaches residual `-228292`. These circuits integrate
the signed carry through the full low `w0-1` output bits. Mixed tables can have
negative constants before arithmetic as well; a `2^21` bias makes them loadable,
and a reversible constant subtraction removes that bias before phase. The
final depth-7 models instead have residual range `[147391,8388606]`, permitting
integration through the 23-bit residual only. An early unverified circuit that
omitted the extended carry was discarded.

## Reversible architectures and lifetimes

### A. Variable-depth local coordinates

[recenter.py](recenter.py) copies only the prefix bits that may be overwritten.
For maximum refinement `K`, it holds `K` copied prefix bits and `K` threshold
flags. Each refinement level changes the upper coordinate bits to the sign
extension of the next centered bit, using one temporary AND. Undoing it costs
one more AND per level; both directions erase temporary products with the
measurement phase correction. No full offset register is materialized.

The component passes 204,800 native cases across bases 5–9 and refinements
0–4, checking the forward coordinate, restored inputs, unchanged flags, clean
scratch, and measurement phase. The final local circuit needs one copied
prefix bit and one threshold flag. The coarse 184-leaf design needs the same
two extra wires but wider polynomial arithmetic. Fewer leaves do not compensate
for that arithmetic.

### B. Shared coordinate with asymmetric prefix refinement

The coordinate stays centered at depth 6, 7, or 8, while individual binary-prefix
leaves refine independently. The lookup recursion stops on equal rows and uses
one AND per non-root internal node. This removes the variable-coordinate wires
but forces all leaves to tolerate the common coordinate's rounding noise.
The depth-6/27-bit-output model demonstrates that refinement alone can fail;
widening derivative precision fixes its numerical accuracy at a resource cost.

Both builders use in-place Horner accumulators, one-bit signed carry counters,
late loading of high constant bits, early HMR deletion of high curvature bits,
deferred allocation of phase headroom, borrowed product/carry workspace, and
replay of curvature to repair its recorded measurement phase. Normalization
flags also supply the parity address and output scaling. The normalized leading
input bit supplies the fixed top bit of the phase word.

### Measured lifetime profile

The tables below use actual builder liveness and the compiled wire-recycling
pass. They do not sum independent component width estimates. Allocator-reserved
high-water marks can exceed the final physical width because a dead wire is
reusable before the Python builder frees its containing array. The compacted
physical circuits and their scratch contracts are what native verification sees.

| Stage entry | Shared live | Local live | Relevant lifetime |
| --- | ---: | ---: | --- |
| Normalization | 27 | 27 | 22 input + 5 normalization flags |
| Initial lookup | 74 | 76 | Residual 23, slope 17, curvature 7; local adds 2 |
| First Horner stage | 86 | 88 | Ten private extension bits + two counters allocated |
| Second Horner stage | 86 | 88 | Three curvature bits freed; three extension bits allocated |
| Phase application | 86 | 88 | 37-bit phase word aliases the normalized top input bit |
| Curvature replay | 86 | 88 | Three extension bits freed before restoring curvature |
| Lookup erasure | 74 | 76 | Arithmetic extension and counters released |
| Denormalization | 27 | 27 | Coefficients and coordinate auxiliaries released |
| **Compiled physical peak** | **90** | **92** | Includes temporary products/carries and lookup recursion |

Two high residual bits are borrowed until the final constant load, then measured
early in reverse. Curvature measurement outcomes are repaired on replay;
delayed residual outcomes join final lookup erasure. Output moves into clean
headroom use AND/CNOT forward and measured-AND deletion backward. Input
normalization retains reversible swaps.

### Component costs of the frozen candidates

| Component | Shared static Toffolis | Local static Toffolis |
| --- | ---: | ---: |
| Normalize / denormalize | 109 / 101 | 109 / 101 |
| Initial coefficient and depth lookup | 267 | 285 |
| Coordinate transform / restore | 0 / 0 | 1 / 1 |
| Curvature Horner / reverse | 105 / 105 | 105 / 105 |
| Linear Horner including delayed lookup / reverse | 417 / 386 | 417 / 386 |
| Affine restore / remove | 8 / 8 | 8 / 8 |
| Output shifts and phase block | 217 | 217 |
| Curvature replay | 33 | 37 |
| Final measured lookup erasure | 37 | 37 |
| **Total** | **1793** | **1817** |

The final local circuit pays 18 lookup gates, 4 replay gates, 2 recentering gates,
and 2 physical qubits beyond the final shared circuit. Its slightly larger
synthesis margin does not repay these costs.

## Approximation versus resource frontier

The full compiled table is [resource_frontier.csv](results/resource_frontier.csv).
The following costs use **static** Toffolis plus the model's exact-margin
synthesis work; they are conservative comparison ceilings, not native scores.
The final native table above uses measured mean Toffolis instead. In these
circuits normalization cleanup saves about 11 mean Toffolis relative to static.

| Construction / schedule | Q | Static Toffolis | Compiled/model cost ceiling |
| --- | ---: | ---: | ---: |
| Final shared, 269 leaves | 90 | 1793 | 183946.52 |
| Shared, constrained to 88 qubits | 88 | 1856 | 185402.82 |
| Shared, constrained to 86 qubits | 86 | 1926 | 187209.12 |
| Final local, 287 leaves | 92 | 1817 | 189774.25 |
| Local, constrained to 90 qubits | 90 | 1883 | 191588.72 |
| Local, constrained to 88 qubits | 88 | 1958 | 193931.19 |
| Shared finer coordinate, narrower coefficients | 88 | 1945 | 192591.91 |
| Shared depth 6, wider derivatives, 217 leaves | 96 | 1921 | 207939.72 |
| Local depth 6–8, 285 leaves | 94 | 2001 | 211195.77 |
| Local 184 leaves, 28-bit output | 100 | 1941 | 219438.95 |
| Local 163 leaves, 29-bit output | 104 | 1942 | 229080.99 |
| Mixed, 1/32 linear tail, streamed cleanup | 106 | 2101 | 249565.29 |
| Mixed, 1/8 linear tail, streamed cleanup | 106 | 2201 | 260165.29 |
| Mixed, 1/2 linear tail, streamed cleanup | 106 | 2846 | 328535.29 |
| Linear, streamed cleanup | 100 | 3121 | 338043.09 |

Schedule grids vary the qubit cap, retained curvature bits, delayed constant
bits, product-bank limit, and deferred phase wires. Feasible and infeasible
attempts are retained in `results/*_schedules.jsonl`; all successful resource
comparisons were compiled and wire-recycled. They include schedules that
prioritize Q over Toffolis. `results/audit.json` counts 362 feasible and 376 infeasible retained schedules.

### Alternative low-workspace cleanup

Two resident ANF monomial banks require
`(2^a-a-1)+(2^b-b-1)` scratch wires. Streaming the larger bank depth first requires
`(2^a-a-1)+(b-1)`, for `a=floor(n/2), b=n-a`, and computes exactly the same number
of temporary ANDs. Every streamed high monomial applies its complete row of
classically conditioned CZ repairs before deletion.

The independent component check exhausts all addresses for random tables at
address widths 3–11: 4,088 native shots pass input, scratch, phase, and equal-gate
checks. Complete streamed circuits also pass 192,512 native linear shots and
119,808 native mixed shots over random inputs and dyadic boundaries, with
repeated measurement seeds. These are shortened checks.

| Same polynomial and arithmetic schedule | Resident cleanup | Streamed cleanup | Static Toffolis |
| --- | ---: | ---: | ---: |
| Linear | 111 Q | 100 Q | 3121 in both |
| Mixed, half linear | 114 Q | 106 Q | 2846 in both |

This preserves a real component improvement even though neither complete
architecture is competitive. The mixed comparison with 119 qubits and 2702
Toffolis uses a different, less constrained arithmetic schedule; it is not the
equal-gate comparison in this table.

## Interpretation and remaining limits

The numerical advantage of a small adaptive table is real: 163–217 leaves can
meet the error requirement at suitable precision. Its arithmetic is more
expensive than the incumbent's, and the final physical peak rises. Narrowing
the coefficients or output instead produces more leaves or a rounding floor.
Mixed degree adds deep address bits and lookup entries while retaining the
quadratic schedule wherever any quadratic leaves remain.

This is a negative result for the **explored** integer-Horner, prefix-lookup,
and register schedules. It is not a lower bound on all adaptive interpolation.
The fitter is heuristic, mixed-degree policies cover three tail sizes rather
than all trees, and no global joint coefficient/scheduling optimizer is claimed.
The 269-leaf shared result also shows why counting nominal table rows overstated
the available interval-selection savings. No superior complete circuit has been
selected, and the incumbent remains the recommended submission.

## Reproduction

Run from the repository root. The environment file disables bytecode writes so
read-only dependencies are not modified.

```bash
source initial_experiments/5/inv-dyadic/env.sh

# All 20 independent fits and exact complete-domain certificates.
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/reproduce.py --models

# Standalone coordinate and streamed-erasure native component checks.
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/reproduce.py --components

# Rebuild both final candidates, assert identical bytes, run shortened checks.
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/reproduce.py --circuits

# Retained wide register-schedule frontiers.
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/search.py \
  initial_experiments/5/inv-dyadic/models/shared7_loose_q27.json --wide
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/search.py \
  initial_experiments/5/inv-dyadic/models/local7_q27.json --wide

# Two unchanged-host verifier runs: 200,540,928 shots each, default tolerance.
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/reproduce.py --full

# Validate saved full native evidence and regenerate CSVs/source hashes.
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/summarize.py
```

Individual models can also be certified with `certify.py models/NAME.json`.
For example, the small mixed construction reproduces with

```bash
"$EXP_PYTHON" initial_experiments/5/inv-dyadic/build.py \
  initial_experiments/5/inv-dyadic/models/mixed_tail1_32_q28.json \
  --limit 106 --keep 0 --deferred 5 --stream-erase
```

The Python executable and compatible C++ runtime paths in `env.sh` are the
existing workspace's environment, not new dependencies installed by this task.
