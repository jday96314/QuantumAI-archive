# Phase-oriented inverse square root

Two fundamentally different complete phase implementations have been built:
a small biased product bank, and direct signed-row phasing that never forms the
completed last product or final polynomial word. Neither improves the incumbent.
The accepted round-4 circuit remains the fallback. All project writes are
confined to this folder. Both candidates pass the complete native workload,
200,540,928 shots each, and independently rebuild byte-for-byte.

Read-only references: round-4 `zero_head_affine_invsqrt.py`, its fixed-top and
affine builder chain, `binary_centered_invsqrt.py`, radix-four multiply-add,
short adder, phase-headroom and ANF cleanup; round-5 `inv-dyadic/RESULTS.md`;
round-4 `PHASE_AND_GROUPED_RESEARCH.md`, Goldschmidt row model and Walsh probe;
the unchanged host inverse verifier.

Baseline: Q=90, mean Toffolis=1782.000090874218, 34 rotations,
cost=182956.5271606425. Its minimum phase-synthesis margin is
1.0839411508124785e-12 turns, distinct from its worst relative-error ratio.

## Exact incremental-phase identity

Let k=21-floor(log2(x)), M=x*2^k, s=floor(k/2), and let z be the
14-bit two's-complement encoding of `(M mod 2^14)-2^13`.
For the incumbent table, define H(A,z,e) as the sum of signed radix-four
rows, each rounded independently, with round-half-up for its positive
magnitude. Then

    B = C1 - H(C2,z,-11)
    P = C0 - H(B,z,-12)
    phase = P*2^s / 2^37.

Writing J=H(B,z,-12), the exact phase factors into C0 and -J contributions.
No whole-product rounding can replace H. Carries between rows need not be
constructed when phasing individual integer rows, because the phase is a
character of integer addition. Conversely, modular wrap in a narrow bank
requires an explicit phase correction unless its modulus maps to whole turns.

The copied coefficient table is unchanged. No segmentation, coefficient fitting,
or general Horner scheduling search was performed. Source:
[model.py](model.py), [components.py](components.py), [build.py](build.py).

## Architecture A: phase the table and a small biased product separately

Use the exact identity

    C0-J = (C0-2^18) + (2^18-J) = A+U.

The independent all-input model proves J in [-261376,259793], so U is in
[2351,523520] and fits an unsigned 19-bit register. A fits 27 bits. This bias
is essential: merely treating a signed 19-bit modular product as an unsigned
phase would introduce a nonintegral phase jump of 2^(19+s-37) on wrap.

The circuit loads A, phases A*2^s/2^37 and measures its 27 table bits. It keeps
the random outcomes for deferred phase repair and frees the whole phase bank.
It then loads C1/C2, computes B, and forms only U. It phases U*2^s/2^37,
reverses the product, restores the curvature needed to reverse B, and repairs
all lookup measurements. There is no completed P register.

The two exponent shifts use 37- and 29-bit banks. The second bank saves eight
phase-word wires, but the first and second phases require 34+29=63 rotations.
Separate table loading also costs 261 extra Toffolis. Releasing all curvature
wires reduces the peak but requires 186 Toffolis to replay them. These costs
are included, rather than treating the table phase as a free operation.

Frozen circuit: [split_q84_k0.kmb](candidates/split_q84_k0.kmb), SHA256
`f65f4e4e6715dd66bf59fb2910b25eb15152a43b54c3a66171abe2df407c25dc`.

## Architecture B: phase signed rounded rows directly

[direct_rows.py](direct_rows.py) never adds the seven final product rows
together. For a Booth digit d_g, its magnitude controls are

    odd = z_(2g) XOR z_(2g-1)
    two = NOT(odd) AND (z_(2g+1) XOR z_(2g))
    sign = z_(2g+1).

The controls odd and two are disjoint. Selected source bits can therefore be
XORed into a temporary row bank without introducing an integer carry between
these alternatives. The independent rounding bit is subsequently added with
an actual reversible increment. Omitting that increment would be incorrect.
For e<0, the row magnitude is exactly

    floor(B/2^(-e)) + bit_(-e-1)(B).

After shifting a nonnegative row word V by s, sign is handled without duplicating
every non-Clifford rotation. For an F-bit word and Boolean sign t,

    integer(V XOR (t*(2^F-1))) - t*(2^F-1) = (1-2*t)*V.

Thus CNOT sign into every phase-bank bit, apply its rotations, undo the CNOTs,
and apply ONE correcting rotation on sign. The coefficient is negative here
because the required product contribution is -J. This identity applies after
exponent routing, so the correcting angle is independent of s.

The seven core row widths are 7,9,11,13,15,17,18. Their private banks each have
ten exponent-headroom wires. The circuit uses 160 bit rotations plus seven
sign corrections for the rows, and 34 rotations for C0: 201 in total. All
rounding increments, row deletion, signed corrections, exponent routing, table
repair and reverse first-stage arithmetic are included.

Each row bank is unshifted, its rounding increment is reversed, and its selected
products are erased by HMR with all appropriate CZ repairs. This clears the
row without reversing a full multiply-add into a separate output accumulator.

Frozen circuit: [rows_q84_k0.kmb](candidates/rows_q84_k0.kmb), SHA256
`eb68f6d40c2467e96b6ee42964226b254c344c1841b3afd5df03ee75a0d7e8fc`.

### Deferred measurement repair

Measuring a computed bit f(a) in the X basis with outcome m introduces
(-1)^(m*f(a)). Measurement alone is not clean uncomputation. C0 outcomes are
retained until final table repair; its packed ANF includes C0 alongside the
original slope and curvature. Curvature outcomes are repaired when the exact
high curvature bits are replayed, before reverse arithmetic. Product erasures
repair every selected control/source pair. These corrections are implemented
coherently or conditioned only on random HMR outcomes, never on measured inputs.
The input normalization is fully reversed.

## Exact error and approximation frontier

[model.py](model.py) is independent of kickmix and the circuit builders. It
uses integer bit length, individual signed row rounding, exact modular 64-bit
phase, and the host reference formulas

    r = isqrt(floor(2^130/x))
    target = ((r+1)//2) mod 2^64
    tolerance = r >> 25.

It exhausts all 3,133,452 original inputs and establishes both phase
decompositions exactly. B is in [32816,130688]. The resulting P is in
[67108897,134217726]. Worst error/tolerance is 0.9796714395489131 at x=958102;
the minimum **absolute** synthesis margin is 19,995,185 units of 2^-64 turns,
or 1.0839411508124785e-12 turns. Evidence: [model.json](results/model.json).

[frontier.py](frontier.py) also evaluates genuinely phase-oriented rounding
changes, with no coefficient refit. Evidence:
[error_frontier.json](results/error_frontier.json).

| Representation | Worst error/tolerance | Model |
| --- | ---: | --- |
| Original independently rounded rows | 0.97967144 | PASS |
| Same quadratic, exact rational products | 0.62666001 | PASS |
| Each final row rounded to nearest 2 units | 1.89478505 | FAIL |
| Each final row rounded to nearest 4 units | 2.93736563 | FAIL |
| Each final row rounded to nearest 8 units | 5.33671443 | FAIL |
| Each final row rounded to nearest 16 units | 9.44606569 | FAIL |
| Each final row rounded to nearest 32 units | 18.98989822 | FAIL |

For the exact rational quadratic, with signed integer coordinate z,

    P_exact = C0 - C1*z/2^12 + C2*z^2/2^23
    phase_units = (C0*2^23 - C1*z*2^11 + C2*z^2)*2^(s+4).

This identity eliminates BOTH intermediate row-rounding errors and passes the
exact all-input check. Its absolute synthesis margin is 295,299,560 units,
1.6008221224300678e-11 turns, about 14.77 times the incumbent margin. It is a
numerical result, not a built circuit or resource improvement. Directly
expanding its bilinear/trilinear bit phases would require many rotations;
the larger numerical margin does not make those rotations free. A compact
transposed or square-oriented implementation remains unresolved.

## Complete resources and live state

The following are compiled physical peaks, not sums of independent estimates.
The mean correction for normalization is about 11 Toffolis. The full host cost
uses the actual native mean, including random measurement corrections.

| Circuit | Q | Static Toffolis | Rotations | Synthesis work T/4 | Task cost |
| --- | ---: | ---: | ---: | ---: | ---: |
| Retained incumbent | 90 | 1793 | 34 | 250.85021 | 182956.52716 |
| Biased product bank | 84 | 2293 | 63 | 472.44564 | 231373.43761, full native |
| Direct row phase | 80 | 2434 | 201 | 1553.87887 | 318150.29471, full native |

The synthesis calculation uses the verifier's accumulated margin exactly:

    eps = pi * minimum_margin_turns / rotation_count
    T_synth = rotation_count * (0.56*log2(1/eps)+5.3)
    score = Q * (mean_Toffolis + (T_base+T_synth)/4).

Both candidates have one base T. Per-rotation eps is 5.405240565581188e-14
for the split bank and 1.694179878764253e-14 for direct rows. The better exact
rational-quadratic margin would save only about 34.3 and 109.3 work units at
these respective rotation counts, even if a circuit could retain the same
gate structure. That hypothetical saving would not close either gap.

| Component, both directions where relevant | Incumbent | Split bank | Direct rows |
| --- | ---: | ---: | ---: |
| Normalize / denormalize | 210 | 210 | 210 |
| Parameter/constant lookup | 267 | 528 | 528 |
| First Horner stage / reverse | 210 | 210 | 210 |
| Last-product arithmetic / reverse or row construction / cleanup | 803 | 750 | 311 |
| Affine restoration / reverse | 16 | 0 | 0 |
| Exponent routing | 217 | 372 | 952 |
| Curvature replay | 33 | 186 | 186 |
| Final measured table repair | 37 | 37 | 37 |
| **Static total** | **1793** | **2293** | **2434** |

Direct row construction does remove 492 arithmetic Toffolis compared with the
incumbent's last stage. Its seven exponent shifts cost 734 Toffolis; C0 costs
another 218. For a c-bit row plus ten headroom bits, the four forward/reverse
shift stages cost (2c-1)+(2c)+(2c+2)+(2c+1)=8c+2 Toffolis. Summing c over the
seven rows gives 734 exactly. Its shifts also add 70 rotations compared with
fixed-exponent row banks. These concrete costs explain why the smaller
arithmetic construction is not a complete resource win.

| Stage | Live at entry, both candidates | Split peak | Direct-row peak |
| --- | ---: | ---: | ---: |
| Normalize | 27 | 27 | 27 |
| Constant lookup | 27 | 72 | 72 |
| Constant phase | 64 | 64 | 64 |
| Derivative lookup | 27 | 62 | 62 |
| First-stage arithmetic | 62 | 68 | 68 |
| Residual / row phase | 45 | 84 | 80 |
| Curvature replay | 45 | 62 | 62 |
| Reverse first-stage arithmetic | 62 | 68 | 68 |
| Final table repair | 51 | 51 | 51 |
| Denormalize | 27 | 41 | 41 |

Machine-readable costs and lifetime events are in `results/*_build.json` and
[resource_frontier.csv](results/resource_frontier.csv). Final-stage narrow-cap
feasibility probes at split/78/keep0 and split/80/keep7 failed the adder's live
workspace check. They are not lower bounds on alternative adders or schedules.

At the measured 84 qubits and 63 rotations, the split circuit would need mean
Toffolis below 1551.36 to reach 170,000, versus 2282.00 observed. Direct rows
at 80 qubits/201 rotations would need below 571.12, while the unchanged stages
outside its row block already cost approximately 1378. Thus merely improving
the tested row arithmetic cannot reach the target. A different exponent-phase
mechanism and a different table schedule would be necessary. These statements
are scoped to the tested implementations, not mathematical lower bounds.

## Verification and frozen evidence

Both complete candidates pass 65,536 boundary/random native model checks with
independent frozen-candidate salts; they include 62,339 distinct inputs.
Evidence: `results/*_check.json`. The exact model covers all inputs, separately.

The row component additionally passes all 11,264 combinations of a 4-bit
source, 6-bit signed coordinate and exponent s=0..10, at four independent
measurement seeds: 45,056 native shots. It checks exact phase, arbitrary signed
rows, all public registers and all scratch. Its full-width un-compacted
standalone component has 71 qubits, 1045 Toffolis and 167 rotations, including
exponent routing and every cleanup. Evidence:
[components.json](results/components.json),
[direct_rows_component.kmb](candidates/direct_rows_component.kmb).

The split-bank candidate **passes all 200,540,928 native shots** under the
unchanged host verifier, salt `inv-phase-split-final-20261024`: 84 qubits,
mean Toffolis 2282.000049575915, score 231373.4376130803.
Evidence: [split_full_verify.json](results/split_full_verify.json).
The direct-row candidate **passes all 200,540,928 native shots**, independent
salt `inv-phase-rows-final-20261024`: 80 qubits, mean Toffolis
2422.999811469906, score 318150.2947084351.
Evidence: [rows_full_verify.json](results/rows_full_verify.json).
Total final experimental native shots: **401,081,856**. Both preserve the same
worst error and absolute synthesis margin as the retained circuit.

Both candidates rebuild byte-for-byte; see [rebuild.json](results/rebuild.json).
[audit.json](results/audit.json) binds candidate hashes to native JSONs, exact
shot counts, source hashes and model worst error. The retained fallback is
[incumbent.kmb](reference/incumbent.kmb), SHA256
`27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1`.
Its copied unchanged 200,540,928-shot evidence is
[incumbent_verify.json](reference/incumbent_verify.json).

## Assessment

**A. Verified improvement:** none; the incumbent remains selected.

**B. Validated research breakthrough:** no complete resource breakthrough is
claimed. Exact incremental phase kernels and deferred cleanup work, but their
complete costs exceed the incumbent.

**C. Scoped negative result:** both independently implemented phase architectures
lose after complete lookup, live-width, signed rounding, exponent routing,
measurement repair and synthesis accounting. Independent row quantization
also fails even at the smallest tested truncation. These are failures of these
representations, not a lower bound on phase-oriented inverse square root.

**D. Unresolved hypothesis:** exact rational evaluation of the existing quadratic
has substantially more numerical margin. A phase-aware middle-product or square
primitive might exploit this without a full high-precision product bank. A
shared exponent-phase mechanism could also change the direct-row result. No
competitive implementation or general resource bound is established for either.

## Reproduction

Run from this directory. The environment disables Python bytecode writes to
read-only references and places temporary files in this experiment. It uses
existing dependencies; none were installed or modified.

```bash
source env.sh
"$EXP_PYTHON" model.py
"$EXP_PYTHON" frontier.py
"$EXP_PYTHON" build.py --mode split --limit 84 --keep 0
"$EXP_PYTHON" build.py --mode rows --limit 84 --keep 0
"$EXP_PYTHON" check_components.py
"$EXP_PYTHON" check.py candidates/split_q84_k0.kmb --salt phase-split-frozen-v2
"$EXP_PYTHON" check.py candidates/rows_q84_k0.kmb --salt phase-row-frozen-v2
"$EXP_PYTHON" reproduce.py

# Two unchanged full host tests, 200,540,928 shots each, two workers per test.
"$EXP_PYTHON" reproduce.py --full
"$EXP_PYTHON" summarize.py
```

The direct host command, if running each separately, is

```bash
"$EXP_PYTHON" "$PHASE_PROJECT/from_host/verifiers/invsqrt/verify_submission_invsqrt.py" \
  candidates/split_q84_k0.kmb --max-workers 2 --salt inv-phase-split-final-20261024 \
  > results/split_full_verify.json
"$EXP_PYTHON" "$PHASE_PROJECT/from_host/verifiers/invsqrt/verify_submission_invsqrt.py" \
  candidates/rows_q84_k0.kmb --max-workers 2 --salt inv-phase-rows-final-20261024 \
  > results/rows_full_verify.json
```
