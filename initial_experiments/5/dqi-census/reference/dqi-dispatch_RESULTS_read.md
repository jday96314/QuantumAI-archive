# Exception-correct compact DQI dispatch

**Reopened under the revised limits. Current status: fixed-schedule classical
model passes all 102,400 host shots; native shared-cofactor work is in progress.
No replacement accepted.**

The new fixed model keeps the incumbent boundaries: initial k=0 slot, block3 at
k=1,4,7,10,13,16,19, then serial k=22..34. Host results: 0 wrong decodes, 403
deferred boundaries, 746 missing-degree history units, 5,329,292 intermediate
checks including 3,178,892 EEA quotient assertions. Metadata is exactly six
completion bits, two six-bit gap tags, and three carry-mode bits. The physical
70-word chart is packed and reconstructed from those bits alone. All-weight
controls (4,608 cases) also pass. Million-case fixed runs are pending.

See [WIDTH_PLAN.md](WIDTH_PLAN.md) for the pre-emission allocation contract:
planned peak 973 Q, with at most six fresh head words. This is an estimate to
validate, not an emitted width. The reused dqi-core files are its checked
trimmed_middle.py and budget_inverse.py only; its full builder is not ACCEPTED.

The revised kill limits are Q>990, projected added T>300k, or >20 wrong host
decodes. The old matrix/Bennett result below is historical and does not terminate
this reopened work.

## Historical matrix/Bennett experiment

The selected round-4 circuit remains read-only and retained, with rechecked SHA256
`4989da5cb633529e8851b9aa0c8c52abeb8990423a3a9463a70f673661a909c0`.
Q=961; mean Toffolis=1,019,865.856484375; 380 phase / 0 value failures in
102,400 native shots; R=1.6968460732519504; cost=1,598,204,229,384.9265.
The combined reference remains 9.3677221531. No archive, verifier, confidence
function, scoring code, chemistry asset, or inverse asset was modified.

## First numeric checkpoint and classical results

Before native resource screening, the autonomous classical model recovered **all
102,400 exact original host cases**, without an exact-Euclid dispatch hint.
Independent million-case streams had **1 and 0 wrong decodes**. These are
classical results, not native decoder acceptance.

| Classical inputs | Wrong decodes | Independent-mask phase failures | Full-state EEA audits |
| --- | ---: | ---: | ---: |
| Exact original host shots, 102,400 | **0** | **0** | 12,094 |
| Law sample, 1,000,000, seed 2026102405 | **1** | **1** | 12,085 |
| Law sample, 1,000,000, seed 2026102406 | **0** | **0** | 12,097 |
| All weights 0–35, 128 cases each, seed 2026102423 | **0** | **0** | 29,698 |

Every case is compared against an independent Berlekamp–Massey locator and its
true error magnitudes. A matching degree-w locator has exactly the w supplied
distinct roots, excluding extra decoded roots. Exact host decoded values imply
zero phase failures under **every** mask, including the original masks. Masks
for the single law-sample failure use a separate RNG; both fixed seeds are
reported. No favorable masks or salts were selected.

The first 1,024 cases of each large run and every weight-control case also
compare every block's Bprev, B, retained remainder, and normalizer against
conventional full polynomial EEA: **65,974 state comparisons**, zero disagreement.
Exact EEA only supplies assertions, never a decoder branch.

The host model encounters 281 boundary deferrals, 23 completions immediately
after a gap, 170 serial gap quotients, and one final-two-slot gap. Nontrivial
patterns include 186 `(1,2)`, 390 `(2,1)`, and 170 isolated `(2)` transfers.
The law streams also exercise `(3)` and multiple gaps.

The sole failure is retained in
`results/law_1000000_2026102405_g3_failures.jsonl`. Independent full EEA diagnosis
is in `results/serial_failure_diagnosis.json`: a degree-three quotient exceeds
the two-slot serial budget. It is below the approximately 20-per-million limit.

## Autonomous local-transfer model and limits

`local_transfer.cpp` retains monic B of degree d, monic Bprev of degree h, band
`C=R[d:N-d]`, and pinv. The incumbent's moment identity generalizes to

```
F = reverse(C)*reverse(B)
    / (1/pinv + z**(d-h+1)*reverse(C)*reverse(Bprev))
```

For g-1 initial zeros, set `beta=F[g-1]` and `W=beta/F[g-1:]`. The reversed
first g+1 coefficients of W give the monic degree-g quotient; `W[g+1:]` is the
next local series. The same recurrence composes a polynomial 2x2 transfer for
all patterns, without inverting a missing pivot. Its determinant is checked
against the product of betas. One shared application and tail update follow:

```
Bprev_new = T00*Bprev + T01*B
B_new     = T10*Bprev + T11*B
pinv_new  = pinv / product(beta)
C_new     = reverse((reverse(C)*reverse(B_new))[2*consumed:] / reverse(B))
```

Blocks allow three degree slots before d=22 and two afterward. Quotients that
straddle the remaining budget are deferred; later blocks use the actual degree
and bank alignment. A zero retained remainder band completes immediately (the
omitted low evaluator polynomial can remain nonzero). At most six
local moment words are read; at most 70 field words store the polynomial state.
Unused exceptional coordinates are implicit zeros.

**This is a variable-degree classical model.** Its metadata and adaptive block
boundaries do not prove a fixed quantum schedule, one-bit physical carry, or
reversible release of exceptional coordinates. Those operations remain
unimplemented. Degree-three gaps also require more alignment information than
the isolated degree-two case. The earlier repair's first-stage hint and
single-gap splice are not used in this decoder.

## Native components

`matrix_transfer.py` emits a uniform reversible multiplication of low coefficient
bank interiors by a unimodular polynomial matrix. A reversible scalar row pivot
handles a zero constant entry; a shared LDU path and high-to-low cross updates
apply every pattern. The monic cubic term is a free XOR.

| Native component | Cases per direction | Total component Q | T per direction |
| --- | ---: | ---: | ---: |
| Common matrix interior, 3 coefficient pairs | 256 | 250 | 1,528 |
| Common matrix interior, 18 coefficient pairs | 4,096 | 610 | 8,638 |
| Six-coordinate to matrix-entry Bennett bridge | 4,096 | 330 | 4,356 |
| Three completion flags into six-bit counter | 4,096 | 17 | 75 |
| Two completion flags into six-bit counter | 4,096 | 16 | 45 |

The 18-pair matrix tests have 1,024 cases each of `(1,1,1)`, `(2,1)`, `(1,2)`,
and `(3)` in forward, inverse, and both composition orders. They match independent
native-field polynomial EEA/convolution, preserve public operands as specified,
leave scratch zero, and repair all random measurement phases.

`matrix_coefficients.py` reconstructs moments from six Padé coordinates, computes
matrix entries, and clears temporary moments. Arbitrary initial destinations
test XOR semantics. Forward/inverse/composition, scratch, and phase checks pass.
Its meaning covers `(1,1,1)` and `(2,1)` only: singular minors, completion,
`(1,2)`, and `(3)` need another head schedule. It is **not a complete block**.

The matrix interior omits grafting consumed moments into the upper coefficients,
bank realignment, completion, and decoder phase scan. Component widths include
their public operands and are **not full decoder widths**.

## Measured resource screen and kill

`profile_transfer.py` compares emitted arithmetic at k=1,4,7,10,13,16,19 with
the incumbent's three cofactor passes, including monomial padding. Counts match
the compiled native components at n=3 and n=18. The empty n=0 interior is pruned.

| Contribution across forward and inverse | Additional T |
| --- | ---: |
| Common interiors versus entire incumbent cofactor passes | 22,536 |
| Two adjacent-minor heads versus recorded/unrecorded regular heads | 31,962 |
| One coefficient generation per block and direction | 60,984 |
| **Partial projection before coefficient cleanup** | **115,482** |
| One further coefficient cleanup per block and direction | 60,984 |
| **Projection including that cleanup** | **176,466** |

Underlying measured head counts are 2,094 T for the six-moment map and 969/936 T
for recorded/unrecorded ordinary continued fractions. The native bridge costs
4,356 T per generation. Partial projected mean T is **1,135,348**; including one
cleanup gives **1,196,332**. These are **estimates assembled from measured
components**, not native complete-circuit measurements or accepted scores.
Upper coefficient grafts, completion, carries, additional coefficient work
across passes, all serial repair, and further singular-pattern heads are omitted.

The bridge has 72 coordinate-input bits plus 258 other live qubits. Even allowing
its coordinates to occupy existing public state, embedding it in 840 state bits
with a six-bit count and three local flags requires **at least 1,107 qubits in
this concrete schedule**, before the usual duplicate block head.

This exceeds the **+90k T kill limit** and Q≤965. Work stopped before emitting a
full decoder or running full native verification. No full decoder controls or
independent native salts are claimed. `results/whole_block_selector_screen.json`
separately prices the naive alternative that emits both complete blocks on every
shot: **+369,600 T**, before controlling their gates, completion, or serial repair.

**These failures concern the explicit matrix/Bennett representation and naive
whole-block selection only.** They do not rule out shared-arithmetic dispatch,
uniform transfers, or a different compact basis.

## Width trim retained for reuse

`completion_counter.py` releases local unary flags under the valid completion
suffix promise. After adding a block's flags into six count bits, flag j equals
`counter >= block_width-j`; recomputing it clears the flag. The inverse regenerates
flags before subtracting. Tests represent every valid incoming count/local suffix
and pass both composition orders, clean scratch, and phase.

This suggests **13 fewer qubits** at the incumbent block-19 peak:
`22 flags -> 6 count + 3 local flags`, or Q≈948 before new transfer workspace.
At scan time, reuse the retained count for rotations and expose the leading monic
B word explicitly: `959 - 35 flags + 12 leading bits = 936`. The original scan
already held a six-bit rotation count. The extra 12 bits avoid expensive binary
degree equality tests at every scan position.

**These are layout estimates, not integrated decoder widths.** The incumbent's
spurious gap flags violate the suffix promise, so this is not a drop-in change.
No full zero test, completion detector, binary-count scan, or reduced-width
decoder was emitted. The measured primitives establish reversible history
compression under the promise, not the remaining integration. They cannot
rescue the measured matrix bridge's width or T cost.

## Reproduction and attribution

Run here, with at most two CPU workers. Exact host cases and probabilities are
read-only host-checked inputs in `../dqi-repair/parameters/`. Seeds were fixed
in advance and both streams reported. The field/gate generator adapts
Apache-2.0 andigu/uncompute revision `64b5a13aa249f1360846657392e6c66817b26d0f`;
source/license remain at `../../3/vendor/uncompute`. Round-4 Codex supplies the
compact identities, six-moment map, field conversion, compiler, and ordinary
arithmetic. Round-5 repair supplies diagnosis and the selector-screen gap block.
Its tiny field tables and the round-4 translator are copied locally for relative
imports; originals are unchanged. All new source and output stays here. No other
concurrent experiment directories were consulted.

```bash
source env.sh
g++ -O3 -std=c++17 local_transfer.cpp -o local_transfer
./local_transfer 102400 2026102405 host 3 > results/host_102400_2026102405_g3.json
./local_transfer 1000000 2026102405 law 3 > results/law_1000000_2026102405_g3.json
./local_transfer 1000000 2026102406 law 3 > results/law_1000000_2026102406_g3.json
./local_transfer 4608 2026102423 controls 3 > results/controls_4608_2026102423_g3.json
"$EXP_PYTHON" matrix_transfer.py --n 18 --shots 4096
"$EXP_PYTHON" matrix_coefficients.py --shots 4096
"$EXP_PYTHON" completion_counter.py --block 3 --shots 4096
"$EXP_PYTHON" completion_counter.py --block 2 --shots 4096
"$EXP_PYTHON" profile_transfer.py
"$EXP_PYTHON" diagnose_failure.py
"$EXP_PYTHON" final_audit.py
```

All files in `candidates/` are named component circuits; no native submission
candidate exists. Acceptance still requires a complete decoder, native zero/
one/weight-35/ordinary controls with independent masks, then all 102,400 default
native shots and independent salts under `flock /tmp/r6-fullverify.lock`.
The resource kill prevents that stage. The objective remains
`R(epsilon_plus) * Q**2 * mean_Toffolis`; no classical failure count or component
measurement was substituted into an accepted score.
