# Round 5 inverse width research

**Stopped under the requested 165k projection criterion. No replacement selected.**
The preserved round-4 circuit remains the only fully verified selection at
**182,956.5271606425**, with 90 qubits. No new circuit projected below 165k.
This is evidence about the schedules implemented here, not a lower bound on
inverse-square-root circuits or on the algorithm families generally.

All new artifacts and writes are in `initial_experiments/5/inv-width`. Native
jobs used at most two concurrent single-worker processes; no GPU was used.
The host verifier, confidence, and scoring were not modified. No new full
200,540,928-shot native run was warranted or launched.

## Baseline reproduction and exact synthesis margin

The preserved baseline is `reference/zero_head_affine.kmb`; its full native
result is `reference/zero_head_affine_verify.json`. The local generator rebuilds
it **byte for byte**, SHA256
`27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1`.
The original round-4 files are untouched.

`results/reproduced_objective.json` independently reproduces the full objective:

| Quantity | Value / evidence |
| --- | --- |
| Q | 90, byte-identical rebuilt circuit |
| Static Toffolis | 1,793 |
| Mean executed Toffolis | 1,782.000090874218, preserved full native result |
| Arbitrary rotations / direct T | 34 / 1 |
| Exact minimum remaining margin | 19,995,185 units of 2^-64 turns at x=2,965,859 |
| Rotation synthesis epsilon | 1.0015592812694554e-13 |
| Rotation and direct T equivalents | 250.85021091069865 |
| Objective | 90 * (1782.000090874218 + 250.85021091069865) = 182956.5271606425 |
| Difference from saved objective | exactly 0.0 |
| Exact worst error/tolerance | 0.9796714395489131 at x=958102 |
| Acceptance evidence | preserved full native PASS, 200,540,928 shots |

The minimum synthesis margin and the worst relative error occur at different
inputs. A sampled worst-relative-error match alone does **not** reproduce the
full objective. New projected objectives below use exhaustive integer margins,
not the larger margins observed in samples.

The baseline additionally passes a new 16,640-shot deterministic domain sample,
including normalization boundaries, interval edges, powers of two, uniformly
sampled inputs, and the known worst input. Its stochastic sample mean differs
slightly from the full mean, as expected. See `results/incumbent_native_sample.json`.

## Actual 90-qubit liveness peak

`profile_incumbent.py` traces allocations and stage boundaries, then applies the
same reset / final-use lifetime rules as wire compaction. The builder's
uncompacted namespace reaches 93 wires; the selected compacted circuit uses 90.
The first 90-wire live peak occurs in the **first radix-four row of the final
Horner product**, inside its low-workspace adder. The trace identifies every
live wire, not just total allocated register sizes.

| Raw wire IDs at that peak | Count | Ownership / role |
| --- | ---: | --- |
| q0..q21 | 22 | in-place normalized input; q21 is also borrowed Booth/carry workspace |
| q22..q26 | 5 | binary normalization/exponent flags |
| q27..q49 | 23 | affine constant remainder; delayed high bits q48,q49 are clean carry workspace here |
| q50..q66 | 17 | updated slope |
| q67..q70 | 4 | retained low curvature bits |
| q71..q73 | 3 | reclaimed high-curvature storage, now headroom / arithmetic scratch |
| q74..q83 | 10 | reserved output headroom, used for products and carries during Horner |
| q84,q85 | 2 | the two Horner carry counters |
| q86 | 1 | Booth magnitude temporary |
| q87..q89 | 3 | extra adder carries |
| **Total** | **90** | |

The nested adder at the first peak uses lhs q79, target
q32..q44 plus q84, incoming carry q78, and borrowed carry workspace
q21,q81,q82,q83,q73,q72,q71,q48,q49. Thus some nominal coefficient/headroom
storage already serves multiple roles. `results/incumbent_liveness.json`
contains all IDs, call stacks, instructions, and arithmetic operands for ten
peak examples. Raw trace line 11273 in the earlier unextended trace was the
first peak; use the recorded line in the current JSON after added state tags.

| Stage | Static Toffolis | New native sample mean |
| --- | ---: | ---: |
| Normalize | 109 | 109 |
| Coefficient load | 267 | 267 |
| Forward Horner, including deferred loads | 522 | 522 |
| Affine restore | 8 | 8 |
| Forward output shift | 113 | 113 |
| Reverse output shift | 104 | 104 |
| Undo affine restore | 8 | 8 |
| Reverse Horner, including curvature replay | 524 | 524 |
| Measured coefficient erasure | 37 | 37 |
| Undo normalization | 101 | 90.010 |
| **Total** | **1,793** | **1,782.010** |

The per-stage sample used the annotated, uncompacted trace and 16,640 shots;
its Q=93 is a profiling artifact, **not** the selected width. All stage counts
are before any rotation-synthesis equivalents. The current normalized totals
are 210 static normalization, 304 load/erase, 1,046 Horner, 217 shift, 16 affine
restoration; the older inherited 94-qubit breakdown does not describe this
selected circuit. Evidence: `results/incumbent_native_stages.json`.

## Incumbent width probes

### Low-ancilla adders, delayed coefficients, and measured curvature release

`probe_workspace.py` and `probe_delay.py` test 423 fixed-table schedules:
limits 70..90, curvature retention 0..4, partial-product caps 4..17 or uncapped,
and delayed constant-bit counts 2..5. The inherited adder selects Cuccaro
prefixes when necessary to meet the requested live-space cap. 71 schedules
build; none reaches fewer than 88 compacted qubits. Invalid cap combinations
are recorded explicitly in the delay probe rather than treated as valid
resource results.

| Schedule | Q | Static Tof. | Projected objective |
| --- | ---: | ---: | ---: |
| Preserved selected | 90 | 1793 | **182956.5272**, full native |
| Best delay=3, keep=4, product=17, limit=88 | 88 | 1842 | 183202.82 |
| Erase more curvature: keep=2, delay=2 | 88 | 1856 | 184434.82 |
| Erase all curvature: keep=0, delay=2 | 88 | 1926 | 190594.82 |

The 88-qubit delay=3 candidate passes 16,640 new native shots with unchanged
worst relative error. Its score is **projected**, using the exact unchanged
synthesis margin and expected 11-Toffoli measurement saving, not accepted from
a short run. Its extra 49 static Toffolis cost more than its two-qubit saving.
The ordinary-adder grid's best result remains the 90-qubit selection.

The existing normalization already shifts the 22-bit input in place. Horner
already updates its loaded coefficient registers in place. There is no
out-of-place mantissa or second full Horner destination to remove in this
implementation. Releasing curvature is a real lifetime change, but replaying
its table costs the gates above. This does not rule out a different in-place
encoding of exponent information or a different multiplication algorithm.

### Segmented phasing and outcome-conditioned carry erasure

`segmented_probe.py` independently implements the requested low/high sequence:
compute and phase low sum; retain carry; uncompute low; compute and phase high;
uncompute high; measure the carry and execute its phase repair **only when the
measurement outcome requires it**. The native simulator counts executed
Toffolis, including the conditional repair.

| Addition kernel | Q including public operands | Static Tof. | Native mean Tof. | Native shots |
| --- | ---: | ---: | ---: | ---: |
| 27-bit whole result | 83 | 158 | 158.000 | 8192 |
| Split at bit 14, measured carry cleanup | 71 | 210 | 195.894 | 8192 |
| 8-bit whole result | 26 | 44 | 44.000 | 8192 |
| Split at bit 4 | 23 | 56 | 51.970 | 8192 |

All restore both arbitrary input operands, clean scratch, and produce the
exact expected phase. A retained *single* addition carry works. The remaining
Horner multiplication carry is not a single bit: exhaustive evaluation of the
selected exact Booth pipeline gives a residual in [-259793,261376].

| Low segment | Carry range | Signed carry bits | Low + carry bits |
| --- | --- | ---: | ---: |
| 10 | [-253,255] | 9 | 19 |
| 12 | [-64,63] | 7 | 19 |
| 14 | [-16,16] | 6 | 20 |
| 16 | [-4,4] | 4 | 20 |
| 18 | [-1,1] | 2 | 20 |

A naive bit-14 split therefore still needs a 20-bit low working word, then a
retained six-bit carry while the low calculation is undone. The updated slope
also remains live. Repeated multiplication / coefficient replay and the
extra variable-shift phase headroom remove the simple one-carry kernel's
advantage for this schedule. No full inverse circuit from this segmentation
was claimed to meet the target. These bounds describe the signed-carry
representation tested here; sign-conditioned carry encodings may differ.
Evidence: `results/segmented_native.json`, `results/segmented_incumbent_carry.json`.

### Dirty borrowed input bits: width does fall, gate count rises

`independent_arithmetic.py` supplies a new controlled MAJ / conditional-UMA
adder. Its unequal-width version initially zero-pads its source. A later probe
replaces padding with arbitrary dirty input bits, then subtracts their restored
high-word contribution. This costs extra additions but needs only one clean
carry bit. Both versions pass 20 native kernel cases / 81,920 shots each,
covering subtraction, incoming rounding carries, random dirty padding,
restored sources and zero residual phase.

Cross-pollinating this into the independent complete linear and serial
circuits gave the following measured changes:

| Circuit | Original Q / Tof. | Dirty-borrow Q / Tof. | Full-margin projected objective after change |
| --- | --- | --- | ---: |
| Bit-serial Horner | 102 / 4159 | 90 / 6119 | 571764.35 |
| Table-heavy linear | 96 / 4466 | 84 / 5544 | 486502.74 |

This probe establishes a working width-reduction mechanism, but its gate
tradeoff is poor for these pipelines. Evidence: `results/dirty_borrow_tests.json`
and `results/dirty_borrow_checkpoints.json`.

## Three independent algorithm families

`family_seeds.py`, `independent_arithmetic.py`, `independent_builders.py`, and
`adaptive_linear_probe.py` do **not** import the incumbent model, coefficient
table, arithmetic, or erasure code. They use fresh Chebyshev seeds, then refit
integer coefficients by coordinate descent and constant-offset minimax repair
against each specified integer rounding pipeline. They use the host's public
kickmix gates and reference/scoring math. The baseline-only carry probe is kept
separate and explicitly marked as incumbent reuse.

The numeric first analytic checkpoint was recorded before native circuit work:
approximately 485k direct, 360k serial, and 475k table-linear. Those initial
analytic schedules were optimistic: unequal-width padding was omitted and the
29-bit Horner/linear rotation estimate undercounted by two. They are retained
as historical estimates, **not measured scores**. The native checkpoints below
supersede them. `results/family_precision_screen.json` records the subsequent
precision reductions and rejected fits.

Each final circuit has all 3,133,452 inputs checked by its exact integer model,
4,096 sampled native phases matched **bit for bit** to that model, and a separate
16,640-shot native domain check of phase tolerance, restored input, and clean
scratch. Full-domain margins use host integer square-root reference values in
units of 2^-64 turns. Projected scores include the actual Q, executed native
Toffolis, direct T gates, and the host rotation-synthesis formula at those
full-domain margins. **None is accepted from these development checks.**

| Independent native design | Q | Mean Tof. | Rotations | Exact max error/tolerance | Full-objective projection |
| --- | ---: | ---: | ---: | ---: | ---: |
| Direct quadratic, separate term lifetimes / phases | 126 | 15393 | 87 | 0.95288749 | 2020351.31 |
| Bit-serial Horner, f=8, widths 27/18/7 | 102 | 4159 | 34 | 0.70576022 | 448079.59 |
| Aligned table linear, f=10, widths 28/18 | 96 | 4466 | 35 | 0.87940900 | 452515.13 |

Detailed evidence: `results/family_final_checkpoints.json`; the five native
circuits, including dirty-borrow variants, remain under `candidates/`.

Concrete obstructions for the implemented designs:

1. **Direct quadratic / segmented term phase.** Keeping the polynomial terms
   in separate lifetimes removes simultaneous c0/slope storage, but the exact
   schoolbook square uses a 28-bit square and a 35-bit curvature product. With
   seven curvature bits and 27 input/flag bits this is already 97 wires before
   adder padding. The native schedule reaches 126, executes 15,393 Toffolis and
   pays for 87 rotations across three phases. Truncating or fusing the square
   could change this result; the present exact-square implementation does not
   approach the target.
2. **Bit-serial Horner.** One carry bit removes a product bank, but serial
   controlled additions execute 4,159 Toffolis before dirty borrowing. A
   27/18/8 fit at f=7 failed exhaustive error tolerance (ratio 1.03081), so this
   candidate uses f=8. The first native circuit also caught a real overflow at
   x=4: a normalized value 134217729 exceeded the 27-bit accumulator. Adding a
   high guard bit fixed the mismatch; the resulting model/native matches are
   exact. Dirty padding cuts 12 qubits while adding 1,960 Toffolis. Streaming
   table coefficient bits would trade these stored bits for repeated lookups;
   that is not a gain with the implemented unary tables.
3. **Table-heavy dyadic linear.** The 10-bit mantissa table has 2,048 parity-
   augmented leaves and costs 2,046 Toffolis just to load. Even idealizing
   Q=70, adding normalization and synthesis nearly exhausts the 165k budget
   before arithmetic, shifting, or erasure. Nonuniform power-of-two-aligned
   intervals help modestly: fine spacing 1/1024 below mantissa 1.875 and coarse
   1/512 above it gives 1,920 leaves and exhaustive ratio 0.93478649. An
   optimistic 70Q schedule, with interval-selection gates free, still projects
   to **287423.98**. Moving the boundary to 1.9375 gives 1,984 leaves and
   **291619.94**. These last two are integer-model/schedule estimates, not
   native candidates. More compressed lookup or a different interpolation
   architecture could change that obstruction.

The shared ideas actually tested across families were conditional carry
repair, separate phase lifetimes, integer-pipeline fitting, low-ancilla adders,
and dirty input borrowing. The aligned interval probe adds the table family's
own independent nonuniform fit. None of these observations implies that the
public leaderboard's stronger solutions cannot be reached by another design.

## Decision and preserved selection

All requested checkpoint categories have been completed, and none of these
designs projects below approximately 165k. Per the stop rule, work stops here
with the **90Q / 182956.5271606425** full-verified selection preserved. The
best new incumbent schedule is 88Q / approximately 183203, not an improvement.
The Q<=70 / <=2000-equivalent target and Q<=60 stretch were not achieved for a
complete inverse circuit. The 71Q segmented *addition kernel* must not be
reported as an inverse-square-root solution.

The overall preserved local reference remains **9.3677221531**. DQI
1,598,204,229,385 and chemistry 708,547,388.9553 are unchanged and outside this
assignment. No archive or submission was replaced.

## Reproduction

Run from this directory; `env.sh` disables bytecode writes and caps math-library
threads. Build the local compactor once, then run development commands:

```bash
source env.sh
g++ -O2 compact_wires.cpp -o compact_wires
"$EXP_PYTHON" profile_incumbent.py
"$EXP_PYTHON" reproduce_objective.py
"$EXP_PYTHON" native_sample.py reference/zero_head_affine.kmb --output results/incumbent_native_sample.json
"$EXP_PYTHON" native_sample.py results/incumbent_tagged.kmx --stages --output results/incumbent_native_stages.json
"$EXP_PYTHON" probe_workspace.py
"$EXP_PYTHON" probe_delay.py
"$EXP_PYTHON" independent_arithmetic.py
"$EXP_PYTHON" dirty_borrow_probe.py
"$EXP_PYTHON" segmented_probe.py
"$EXP_PYTHON" family_seeds.py
"$EXP_PYTHON" screen_families.py
"$EXP_PYTHON" independent_builders.py --family serial
"$EXP_PYTHON" independent_builders.py --family linear
"$EXP_PYTHON" independent_builders.py --family direct
"$EXP_PYTHON" finalize_checkpoints.py
"$EXP_PYTHON" independent_builders.py --family serial --dirty
"$EXP_PYTHON" independent_builders.py --family linear --dirty
"$EXP_PYTHON" adaptive_linear_probe.py
```

`results/dirty_borrow_checkpoints.json` applies the original family's exhaustive
synthesis equivalent to the new measured dirty-borrow gate count; their integer
pipelines are unchanged. Sample scores printed by `native_sample.py` use only
sample margins and are diagnostic; use `family_final_checkpoints.json` for the
complete objective projections.

The following full-default verifier command is documented for future candidates
only. It was **not** run for the losing probes:

```bash
flock /tmp/r6-fullverify.lock "$EXP_PYTHON" \
  ../../../from_host/verifiers/invsqrt/verify_submission_invsqrt.py \
  candidates/CANDIDATE.kmb --max-workers 2
```

## Reuse and licenses

The local incumbent generator, short/Booth adders, monomial erasure, workspace
search helper, and C++ wire compactor are copied from the read-only round-4
Codex artifacts, with module docstrings preserved. The local env also makes
prior round-3 / round-2 utility imports and the host kickmix package available
read only. The host kickmix license remains at `from_host/kickmix/LICENSE.txt`.
No `andigu/uncompute` source is imported or copied for these inverse circuits.
The prior DQI decoder's Apache-2.0 attribution is retained in its untouched
artifacts; its upstream is https://github.com/andigu/uncompute at revision
64b5a13aa249f1360846657392e6c66817b26d0f. New independent arithmetic and family
sources are identified above and contain no copied incumbent coefficient code.
