# Round 5 / dqi-phase — completed component experiment

**Decision: preserve the round-4 Codex-selected decoder. No replacement is
promoted and no new full-task cost win is claimed.** The 11-Toffoli `F` point
body works, including arbitrary per-shot classical masks and root gating.
The implemented reversible preparation schedules lose under `R * Q^2 * T`.
The direct half-degree Psi transplant loses even when its missing evaluator
conversion is assigned zero cost in a favorable rejection screen.

All new files are inside this directory. Prior rounds and host files are
read-only. No locator recurrence or exception repair was implemented or
changed. The existing round-4 classical cofactor model is called only as a
read-only control for the phase identity.

## Preserved reference and exact attribution

Selected artifact:
`../../4/codex/candidates/cubic_trace_product_n70_m4095.kmb`.
SHA256: `4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0`.
The binary digest was independently reread and matched to the existing full
native report. The combined reference remains **9.367722153134604**, with
inverse **182956.5271606425**, DQI **1598204229384.9265**, chemistry
**708547388.9553173**.

The existing DQI native result is **102400 shots, 0 value errors, 380 phase
errors**, Q=961, mean T=1019865.856484375. Re-evaluating the host's exact
Clopper–Pearson definition gives epsilon_plus=0.006777580874382957 and
R=1.6968460732519504. That is existing round-4 evidence, not a new validation.

`baseline_audit.py` reads the selected tape and scan wrapper, without rebuilding
the locator:

| Part | Exact static Toffolis |
| --- | ---: |
| Compact transform, scaling and its inverse | 649802 |
| Cofactor unpadding and restoration | 9792 |
| 4095 point bodies at 88 | 360360 |
| Total | 1019954 |

The lower native mean comes from runtime classical zero-mask skips. The tape
contains 68957627 core operations. The current kernel has 36 quantum input bits
(root, prior value, derivative), 42 scratch bits, Q=78 locally: a 42T field
product plus a 46T gated cubic-trace inverse phase. The complete circuit's
961Q peak is the authoritative width, rather than a sum of component peaks.
See `results/baseline_audit.json`.

## Independent phase and mask oracle

The host field is GF(4096), polynomial **0x1053**, primitive element **2**,
little-endian coordinates. The disclosed/compact arithmetic field uses
**0x10eb**. The existing reversible isomorphism was independently checked on
all 4096 field elements and on 4096 random products using a separately written
carryless multiplier and polynomial reduction.

Position j corresponds to **beta_j = 2^(-j)**. With the standard evaluator
and constant-normalized locator,

```
S_i = sum_j e_j * 2^((i+1)j), i=0..69
sigma(z) = product_j (1 + 2^j z)
omega = sigma * S mod z^70
Psi * sigma_prime = 1 mod sigma
F = omega * Psi mod sigma
e_j = omega(beta_j)/sigma_prime(beta_j) = F(beta_j)
```

There is no extra position factor in this standard evaluator convention.
The **selected compact state is a different representation**. After unpadding
it has monic `B` and `A = pinv * Bprev`, with

```
e_j = beta_j^70 / (A(beta_j) * Bprime(beta_j)).
omega_monic = x^70 * inverse(A) mod B
F = x^70 * inverse(A) * inverse(Bprime) mod B.
```

`control_identity.py` calls the unchanged round-4 state model and checks this
against syndromes constructed directly from known errors. Generic degree gaps
are explicitly reported; they are not repaired or silently counted as passes.
`oracle.py` does not use a locator recurrence: it constructs locators from
known errors and checks evaluator/Psi/F identities with independent polynomial
arithmetic. It passes 216 split-locator cases covering every degree 0–35,
3780 root checks, 884520 mask words, and 35 general squarefree polynomial cases.
It covers zero masks, all-one masks, degree zero, a squarefree zero root, and
repeated-root rejection. Repeated-root reciprocal derivatives do not exist;
no extension is claimed for their mathematical inverse.

Every coordinate mask M has a unique trace-dual mu satisfying
`popcount(M & e) mod 2 = Tr(mu*e)`. The native coordinate-to-dual columns are
`[2089,1044,522,261,130,65,32,16,8,4,2,1]`.

Choose U(x)=sum_(i=1)^70 u_i*x^i interpolating mu at positions 0..69. Then

```
residual_mu[j] = mu[j] + U(2^j)
phase = sum_i Tr(u_i*S_i) + sum_errors Tr(residual_mu[j]*e_j).
```

This is valid for arbitrary per-shot classical masks. **Exactly 70 skips are
guaranteed.** The other 4025 independent uniform masks each have probability
1/4096 of zero residual, so the expected total is 70.982666015625 skips and
4024.017333984375 active points. Gu's fixed-mask 71st skip cannot be hard-coded.

`mask_scan.py` emits classical interpolation, the syndrome-linear Z phase,
residual mask conversion, root-gated scan, and the full reverse classical
mask transformation. It uses 840 extra **classical bits**, zero extra quantum
qubits for interpolation, and no Toffolis for mask preparation. At 4095
positions the implementation emits 41942204 interpolation/residual CX
instructions, 524160 basis-map CX/SWAP instructions, 136850 residual-coordinate
map CX/SWAP instructions, and 1610 classically controlled syndrome Z instructions.
These operations are included; the official objective assigns them zero
Toffoli/quantum-width cost. This is a substantial instruction overhead even
though it is uncharged by this objective.

Both native prepared-polynomial scans pass 512 shots at M=128 across degrees
0–35, with all **70 syndromes**, polynomial sources, original classical masks,
quantum scratch and phase checked. Static scan costs are 638T for F and 1334T
for Psi. Executed Toffolis are independently checked against each shot's
nonzero residual masks. These circuits have supplied polynomial registers;
their large public width is not a full compact-decoder width claim.

## Point kernels and near-term controls

Each point body passes 61440 native component shots: all 4096 first-operand
values against the 12 basis masks, all 4096 masks, nonroots and zero operands.
Inputs, classical masks, scratch and exact phase are checked. These are
component tests, not shortened full-decoder acceptance tests.

| Body | Toffolis | Local Q / public Q | Result |
| --- | ---: | ---: | --- |
| Selected inverse-product, cubic basis | 88 | 78 / 36 | Fully verified incumbent body |
| Division-based Forney | 80 | 80 / 36 | Requires omega, absent from compact layout |
| Earlier quadratic-tower inverse product | 92 | 80 / 36 | Loses |
| Prior joint cubic cover | 87 | 78 / 36 | Banked component only |
| Prepared omega/Psi bilinear pairing | 23 | 48 / 36 | 11 root-test + 12 gated bilinear terms |
| Prepared F | 11 | 36 / 24 | 11 root-test; mask parity and conditional CZ are Clifford |

The 80T body cannot consume the selected cofactor as an evaluator. A direct
per-point adapter would need a complete 102T inverse lifetime to obtain the
numerator, giving 182T/point; even at unchanged width its optimistic cost is
1.37734 times the incumbent. With a magically available evaluator, 80T alone
would improve total T by only 3.2%, before conversion.

The 87T point body is reproduced from the earlier frozen cover, with source
and parameter attribution. Its existing **full** native result was 102400
shots, Q961, mean T1015772.10875, 406 phase errors, R1.7223319602461098 and
cost1615697003965.9805, **1.0945% worse** than the selected cost. That difference
is finite-sample R evidence, not proof that its identical phase function is
less accurate. No hash/salt tuning or equivalent circuit selection was done.
No new full candidate is retained from these controls.

## Reversible Psi workspace and lifetime

`psi_schedule.py` adapts the disclosed half-degree pass to an explicitly live
compact prior-cofactor bank. It rewrites sigma's coefficients as square roots,
so sigma=E^2+z*O^2, pads both halves by x^m with m=floor((35-L)/2), runs the
half-degree Euclid tape, seeds a Bezout cofactor, and reverses the Euclid tape
interleaved with cofactor updates. With AE+BO=1, Psi=B^2+z*A^2. All phases are
then undone in the opposite order, including padding and square roots.
The compact completion counter supplies m by a wiring shift; the standalone
component takes those counter bits as inputs. Materializing the virtual
monic coefficient uses CNOTs controlled by the completion flags and is reversed
before the compact core's inverse. It is not an additional field product.

**Measured:** 44465T forward, **88930T build and erase**, **1383 live qubits**.
The native KMB width equals the allocator peak. Exact replay-ID holds are
recorded in `results/psi_liveness.json`; this is not inferred from polynomial
dimensions alone.

| Co-live ownership at the peak | Q |
| --- | ---: |
| 70 existing compact field words, including the scale word | 840 |
| Materialized virtual locator word | 12 |
| 35 completion flags + 6 counter bits | 41 |
| Emitter parity pivots | 3 |
| Six-bit exceptional-event record + pre-swap flag | 7 |
| New 36-word Bezout bank, including EXT | 432 |
| Exact journal scratch replay pool | 48 |
| **Peak** | **1383** |

Forward T attribution: padding1416, pre-swap227, Euclid F15201,
seed390, reverse-with-cofactors R27231; setup is Clifford. Each has its
corresponding inverse. All 35 prior-cofactor words remain live. They cannot
serve as Gu's zeroed syndrome tail: they contain the magnitude information.
An initial attempt to allocate the new cofactor bank without reserving the
journal's physical scratch IDs was rejected by the allocator. Its log is
retained; the corrected schedule holds those IDs before allocation and
successfully replays them.

Native forward reciprocal and complete roundtrip tests pass **4096 shots**,
all degrees0–35, random nonzero locator scales, random retained prior words and
random retained flags. A separate 256-shot forced-exception diagnostic
(zero polynomial, constants, repeated roots, zero derivative, arbitrary monic
polynomials, squarefree zero root) found no roundtrip value, scratch or phase
errors. This is restoration evidence, not a claim that Psi exists on nonsquarefree
inputs or that all rare Euclidean classes are solved. The inherited Euclidean
support windows/limited degree-defect handling remain unchanged.

The schedule restores both cofactor banks exactly on tested cases; conjugating
it by the unchanged compact transform would consequently restore its 70
syndromes. That full composition was **not** integrated or native-accepted.

## Reversible F transforms: implementation and failed budgets

An in-place F map is not free because multiplication by Psi is mathematically
invertible. `modular_schedule.py` implements a complete polynomial XOR MAC:
compute the convolution, reversibly divide by the monic quantum modulus,
copy the remainder, reverse division, and erase the convolution. The in-place
omega-to-F construction computes F into a word bank, swaps banks, and erases
old omega using `omega_old = F*sigma_prime mod sigma`; its inverse repeats the
opposite operations. The old-omega bank is retained until the inverse even
when the identity does not hold, so exceptional inputs still restore.
Derivative construction and erasure are explicit Clifford copies.

The schoolbook control costs 202860T per degree35 modular XOR MAC and 811440T
for the F build/erase lifetime, Q2931. A pointwise split-algebra Fermat inverse
was also implemented as compute/copy/uncompute, with all 21 intermediate
polynomials accounted for. It costs 6256320T at degree35 and has a 9651Q
scratch requirement in the four-array resource harness. No claim that
polynomial squaring modulo a quantum modulus is Clifford was made.

To avoid a straw-man rejection, `fast_modular.py` replaces schoolbook work with
Gu's evaluation/interpolation products and round-4 reversible triangular
middle products. A degree35 modular MAC costs **23352T**, Q2511 including
1680 public input/output bits. This is **native checked at degree35**, with
arbitrary quantum monic moduli and arbitrary incoming target words, including
all source, scratch and phase restoration. Degrees2/3 are independently
cross-checked too. The degree35 gate attribution is

```
2 * 42 * (69 convolution products + 141 triangular-division products
          + 68 low-remainder products) = 23352.
```

The complete in-place F lifetime costs **93408T**, with an exact emitted
allocation peak **2931Q**. The two-MAC-per-direction choreography is also
native tested at degree3 on arbitrary input data. Its peak comprises
1260Q for sigma/omega/Psi, 420Q retained old-omega workspace, 420Q derivative,
828Q product/quotient workspace and 3 emitter pivots. The allocation/gate
count at degree35 is an emitter profile; it is not mislabeled as a full
native DQI result. Out-of-place F costs 46704T for build and erase, but retains
another 420Q F bank and still peaks at 2511Q in this arithmetic component.

Further costs remain before either is a compact decoder: conversion
`A -> x^70/A mod B`, variable-degree routing or padding, and retaining/erasing
all conversion workspace. Psi from the half-degree representation can have a
coefficient of degree35. For L35, reducing that extra coefficient while
retaining its reversible quotient costs another 35 field products in each
direction (2940T total). For smaller L there is a degree-dependent adapter.
Those costs are **not assumed free in any accepted implementation**.
The budget below omits them only to reject already-losing routes favorably.

## End-to-end objective screen and obstruction

`budget.py` uses the actual host CP routine and the preserved R. Static
non-point work is 659594T; uniform classical residual masks give the expected
active-point count above. These are optimistic analytical screens, not
assembled-circuit or new verifier scores:

| Screen | Q used | Mean T estimate | Cost / incumbent |
| --- | ---: | ---: | ---: |
| Psi, full build/erase, zero-cost evaluator conversion | 1383 | 841076.3987 | **1.7080083** |
| F, fast in-place omega transform | at least2931 | at least886196.1907 | **at least8.0829834** |
| F, fast out-of-place omega transform | at least2511 | at least839492.1907 | **at least5.6197922** |

F rows use the arithmetic component's width alone, so the full decoder's
flags, retained degree-35 quotient and missing adapters can only worsen these
figures. They are fixed-degree35 component screens, not distribution-weighted
native whole-circuit measurements. The direct Psi screen alone is decisive:
it would require Q<=1058.223 just to break even at the preserved R, and
Q<=946.503 to achieve1.25x. Its implemented peak is1383. Even the unphysical
R=1 limit gives1608717575917.9885, above the incumbent1598204229384.9265.
With zero failures in102400 shots, host CP would instead give an optimistic
Psi cost1916401647238.6265.

At **Q961**, achieving1.25x permits at most815892.6852 mean Toffolis. A 23T
scan allows **63746.2865T for every preparation and inverse**, less than the
88930T Psi lifetime alone. An 11T scan allows112034.4945T for all preparation,
so measured Psi leaves only **23104.4945T** for F construction, inversion,
all compact-layout conversion and their reversals. The implemented fast F
lifetime is93408T before those adapters. At Q1000 the total preparation
allowance for the11T scan falls to49635.8378T, already below measured Psi.
This is an arithmetic break-even calculation, not an arbitrary qubit cap.

A useful future target is a low-workspace reversible map directly on the
existing prior bank,
`A -> x^70/(A*Bprime) mod B`, retaining B and restoring all core flags and
syndromes. On a split squarefree algebra this is a pointwise involution
(zero inverse totalized to zero), but that fact supplies neither a reversible
circuit nor a free erasure. No low-width implementation within the needed
~112kT total preparation allowance was found here. The current working
whole-task choice remains the selected88T inverse-product decoder;87T is a
reproduced component, while11T and23T require the charged preparation.
This is an obstruction for the implemented schedules, not a lower bound
against every possible reversible polynomial algorithm.

## Reproduction and artifacts

From this directory:

```bash
source env.sh
"$EXP_PYTHON" baseline_audit.py
"$EXP_PYTHON" oracle.py --trials 6
"$EXP_PYTHON" control_identity.py
"$EXP_PYTHON" point_kernels.py
"$EXP_PYTHON" psi_schedule.py --shots 4096
"$EXP_PYTHON" exception_audit.py
"$EXP_PYTHON" modular_schedule.py
"$EXP_PYTHON" fast_modular.py
"$EXP_PYTHON" mask_scan.py --m 128 --shots 512
"$EXP_PYTHON" budget.py
"$EXP_PYTHON" audit.py
```

Builders write only local `candidates/`, `parameters/` and `results/`.
`psi_schedule.py --reuse` and `mask_scan.py --reuse` rerun native simulation
using local KMB files. `fast_modular.py --reuse-native-mac` reuses the saved
unchanged MAC checks, tests F choreography, and recomputes emitter profiles.
`candidates/` contains the point kernels, Psi forward/roundtrip, native mask
scans and modular-transform components. `results/` contains their native
component JSONs, exact emitter budgets and liveness records.

No full candidate meets the requested15–20% integration screen, so no new
102400-shot full decoder run or independent-salt final acceptance was attempted.
A future candidate must still pass the full native102400 shots and independent
validation salts; these component tests do not waive that requirement.
No inverse or chemistry candidate, verifier, confidence setting or scoring
rule was changed.

## Attribution

The external disclosed architecture is
<https://github.com/andigu/uncompute/blob/main/tasks/dqi-reed-solomon/circuits/discovered/README.md>.
The implementation used is the existing read-only Apache-2.0 vendored revision
`64b5a13aa249f1360846657392e6c66817b26d0f`. Its README and license are retained
under `reference/`. `psi_schedule.py` adapts its half-degree Euclid/cofactor
choreography; `fast_modular.py` uses its field-product emitter and round-4
middle-product helpers. Other kernels are read-only round-4/earlier local
components with their original attribution. The independent host-basis oracle,
runtime mask schedule, liveness integration and objective screens are local
experiment code. `results/audit.json` and `parameters/source_manifest.json`
record source and artifact hashes.
