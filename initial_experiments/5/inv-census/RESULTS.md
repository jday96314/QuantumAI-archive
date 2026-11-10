# Inverse architecture census

**No replacement found. The requested1.5× improvement is not achieved.**
The selected, full-verified circuit remains Q 90, mean Toffolis1782.000090874218,
34 rotations, cost 182956.5271606425. It rebuilds byte for byte locally. No new
full verifier run or selection/archive change was made.

Phase 1 compared **eight distinct architectures before building**, including
peak register tables, Toffoli estimates, rotations, synthesis margins, objective
ranges and key uncertainties: [PHASE1.md](PHASE1.md). Phase 2 built both top
exponent mechanisms and a further circular-rotation refinement. Each losing
implementation was examined beyond its first draft. Phase 3's qualified and
conditional build decisions are in [BUILD_PLAN.md](BUILD_PLAN.md).

All writes stayed in this directory. Numerical-library threads were capped at 1;
no more than 2 single-thread jobs ran concurrently, below the CPU cap 3. No GPU
or dependency installation was used. Host code, other agents' files and prior
artifacts were read-only. The supplied combined reference 9.3677, inverse 182956.5272,
DQI 1598204229385 and banked chemistry 708156073.7876 are preserved.

## Numeric checkpoint and complete compiled resources

The first checkpoint established all 3,133,452 integer inputs and 15360 native
component shots per top mechanism. Borrowed headroom was 358 Toffolis/34 rotations;
exact exponent-angle multiplexing was 0 Toffolis/307 rotations. Both preserved the
incumbent's exact phase and synthesis margin. Complete lifetime compilation then
showed that arithmetic scratch prevents the anticipated width saving.

All numbers below are **development projections**, except the preserved incumbent.
MeanToffolis are measured in 52736 native shots per complete new circuit, with two
independent measurement seeds. Synthesis uses the **exhaustive** margin, not the
larger margin observed in native samples.

| Complete circuit | Q | Static Toffolis | Sample mean Toffolis | Rotations | Synthesis work T/4 | Projected objective |
|---|---:|---:|---:|---:|---:|---:|
| Preserved incumbent, full verified |90|1793|1782.000091 full|34|250.850211|**182956.527161 full**|
| Input-borrowed conditional headroom |90|1934|1922.978952|34|250.850211|195644.624640|
| Circular27-bit result with wrap correction |90|1854|1842.978952|37|273.343905|190469.057121|
| Exponent-Walsh rotation selector |88|1639|1628.026927|307|2399.218176|354397.569042|

Evidence: [audit.json](results/audit.json), [complete_model_matches.json](results/complete_model_matches.json).
All three frozen candidate files rebuild byte for byte. They are retained under
`candidates/*_best.kmb` as research artifacts, **not accepted replacements**.
The target is 121971.0181; the1.3× build-plan cutoff is 140735.7901.

Each architecture received 63 new lifetime/arithmetic schedules (caps72,76,80,84,
88,90,94; retained curvature0,2,4; product-bank limits8,12,17). Each has 24 feasible
and 39 infeasible schedules; minimum compiled width 88. There are 189 attempts total.
These vary a new headroom-free lifetime, rather than repeating the old 423 incumbent
schedules. Historical path fields in grid rows name transient build destinations;
only the final `*_best.json` and audit identify retained complete candidates.

## A — conditionally clean headroom: valid, but arithmetic still sets the peak

For normalized M, the low k input bits are zero and s=floor(k/2). Append input
bits9..0, reversed, to the 27-bit polynomial to make a 37-position cyclic word.
The last s positions are zero, so rotating left by s inserts exactly the required
zeros. The first s head positions now contain polynomial bits. Phase those only
when s>j; omit the remaining positions because they contain dirty input bits.
Undo the permutation before reversing arithmetic or erasing tables.

The independent integer model proves the complete masked phase and inverse
rotation on every original x. The kernel passes every reduced4-bit output/6-bit
input instance satisfying the zero-bit promise, at two seeds, plus 12267 distinct
full-width polynomial inputs. Frozen evidence has3072 reduced native shots and
12288 full-width shots. The complete circuit additionally matches the integer
phase **bit for bit** on 26148 distinct original inputs at two seeds, 52736 shots.
Restored public input and zero scratch are checked on every shot.

The phase kernel itself costs 358 Toffolis instead of 217:288 routing swaps plus
70 comparison/product gates in this implementation. The complete best draft's
static count is exactly1793-217+358=1934. Removing all permanent headroom does
not remove the row-product and carry storage it previously supplied. Caps below 88
fail inside the arithmetic; e.g. an80-qubit cap already has 80 live wires before
the nested adder requests carries. Erasing more curvature costs replay gates.
This resolves the proposed lifetime's key uncertainty against a resource win.
It does not rule out other conditionally clean arithmetic primitives.

## A refinement — circular result without borrowing input

Rotate the 27-bit result cyclically by s. A wrapped bit at position j<s needs
weight2^(27+j) instead of 2^j. Add a controlled phase of
`2^j*(2^27-1)/2^37` whenever s>j, and reverse the rotation.
This removes the borrowed-headroom routing and lowers the kernel to 278 Toffolis,
with 37 rotations. Its exact identity is certified on the full integer domain,
and the same reduced/full native component and complete-circuit checks pass.
The resulting complete circuit still peaks at 90 and projects above the incumbent.
This is an implemented second draft of the no-headroom idea, not an untested
suggestion to fix the first draft.

## B — exponent rotations: sparse-selector uncertainty resolved

For each polynomial bit, express2^s on s=0..10 in the Walsh characters of the
four exponent bits. Unreachable s=11..15 are free extension values. The best
real-angle extension needs 11 nonzero characters. A CNOT-only identity
`p*chi(t)=(p XOR t)-t` implements the selector, combining exponent-only corrections
across the27 result bits. Actual circuit count is 307 arbitrary rotations.

[exact_walsh_support.json](results/exact_walsh_support.json) checks all 4368
11-column bases exactly:1024 are nonsingular and all have support 11. A further
[modular certificate](results/modular_walsh_certificate.json) checks all 8008
10-column sets. For every rank10 set, its primitive integer left-null vector v
has v dot(1,2,...,1024) with2-adic valuation at most6. Thus allowing whole-turn
angle lifts cannot represent a result-bit selector using10 or fewer characters
when its phase denominator is2^n,n>=7. Our bit denominators have n=11..37.
Integer matrix products certify the reconstructed null vectors; modular rank
uses a prime larger than every possible minor, avoiding floating rank claims.

This rules out the <=35-rotation hope for this separate-bit exact Walsh
architecture, including whole-turn lifts. It is not a lower bound on general
phase circuits, approximate multiplexing or alternative joint bit encodings.
Its native kernel and complete circuits pass the same development checks as A.
The full synthesis cost2399.218176 explains its loss despite lower Toffolis.

## Additional numerical result — narrow fused quadratic correction

Rather than stop at the losing exponent drafts, the arithmetic alternative was
screened in 30 exhaustive integer pipelines, each over3,133,452 inputs. Keep the
linear term exact, but round the coordinate only inside the quadratic term.
With z=(M mod2^14)-8192, u=floor((z+16)/32), and h=round(u^2/64), evaluate

```
num = c0*2^23 - c1*z*2^11 + c2*h*2^16
phase_units = num * 2^(s+4) modulo2^64
```

This pipeline **passes** with worst error/tolerance0.8046964363863656 at x2424816,
and minimum margin 137902011 units=7.475682995815996e-12 turns. The unsigned
absolute coordinate needs 9 bits, including endpoint256; its rounded square h
needs 11 bits. A signed 9-bit u would overflow: it needs a tenth sign/guard bit.
The wider-accuracy q=4/square-round16 variant passes at ratio 0.6964637323 and
margin 237162772 units. Coordinate rounding by six bits fails at ratio 1.1173053262.
No constant/table refit was needed for these checks.

[The frontier](results/curvature_frontier.json) is numerical evidence only.
No fused-square native circuit, physical width, gate count or objective is
claimed. [BUILD_PLAN.md](BUILD_PLAN.md) specifies the precise kernels, carry and
rounding checks, coefficient liveness problem, and gates needed to qualify it.
A hypothetical Q 80,T 1400,R 50 implementation would cost 140233.8434, barely1.3×;
keeping34 rotations would cost 130987.1389. Those are conditional budgets, not
compiled projections. Table replay or keeping the slope resident can erase them.
This architecture remains unresolved and has not been rejected.

## In-place exponent result

The map E=M for k=0, otherwise E=M-2^21+2^(k-1), is injective over the entire
original domain. When E's top bit is zero, k=ctz(E)+1 recovers the exponent.
The all-input model verifies packing, decoding and restoring M. It is not yet a
native normalization design: the low marker changes the local coordinate, and
on-demand decode/repair must be costed. No five-qubit physical saving is claimed.

## Verification, selection and limits

- Each frozen complete candidate passes 52736 exact native model-match shots,
  using two fixed independent measurement seeds, on 26148 distinct x. This is
  a development check, not the host's full 200540928-shot acceptance workload.
- Each component has 1360 exhaustive reduced cases at two seeds, 3072 padded
  native shots;12267 full-width distinct cases give 12288 further shots.
- Separate 16640-shot host-domain diagnostic samples pass for all three complete
  circuits, including normalization and lookup boundaries and the known worst x.
- The exact all-input numerical margin is19995185 units at x2965859 for all three
  exponent designs. Worst relative error0.9796714395489131 occurs at x958102.
- No salt or hash search was performed. Measurement seeds are fixed diagnostic
  labels. No new full verifier run was justified by a competitive candidate.
- The copied incumbent's preserved full PASS is200540928 shots. Its SHA256 and
  byte-identical local rebuild are both
  `27b14e544481db8eae7bf7f7e61bef93660ada06585cb5e398d3a18af9b50af1`.

## Reproduction

Run from this directory. All commands write only locally. No installs required.

```bash
source env.sh
g++ -O2 compact_wires.cpp -o compact_wires
"$EXP_PYTHON" prototype.py
"$EXP_PYTHON" exact_walsh_support.py
"$EXP_PYTHON" modular_walsh_certificate.py
"$EXP_PYTHON" build_experiment.py --mode borrowed --grid
"$EXP_PYTHON" build_experiment.py --mode circular --grid
"$EXP_PYTHON" build_experiment.py --mode walsh --grid
"$EXP_PYTHON" check_complete.py
"$EXP_PYTHON" curvature_frontier.py
"$EXP_PYTHON" audit.py
"$EXP_PYTHON" native_sample.py candidates/borrowed_best.kmb --output results/borrowed_first_sample.json
"$EXP_PYTHON" native_sample.py candidates/circular_best.kmb --output results/circular_native_sample.json
"$EXP_PYTHON" native_sample.py candidates/walsh_best.kmb --output results/walsh_native_sample.json
```

The full-default command below is documented for a future competitive candidate;
it was **not run** on these losing probes:

```bash
flock /tmp/r6-fullverify.lock "$EXP_PYTHON" \
  ../../../from_host/verifiers/invsqrt/verify_submission_invsqrt.py \
  candidates/FUTURE_CANDIDATE.kmb --max-workers 3
```

`modify_local.py` records the one-time checked transformation of copied baseline
sources; it is provenance, not an idempotent reproduction command. Rebuilds use
the already transformed local sources. Grid outputs retain infeasible attempts.

## Sources and reuse

Gu et al., *Discovering state-of-the-art quantum circuits with large language
models*, Supplement III.1–III.3, https://arxiv.org/html/2610.22367v1, informed the
aligned intervals, cyclic shared coordinate, narrow direct accumulation and
measured table erasure comparison. Their absolute-error task resources were
not transplanted to this stricter relative-error task. No upstream
`andigu/uncompute` source was copied or imported.

The local polynomial/affine builders, Booth accumulator, short adder, monomial
phase-repair erasure, wire compactor and native diagnostic are reused from the
read-only inv-width artifacts (which attribute earlier round 4/round 3/round 2
sources). Their module docstrings are preserved; modifications are local.
The new exponent kernels, exact support certificates, independent integer
models, complete model checker and audit are census work. Kickmix and verifier
APIs are read-only dependencies; the upstream Apache-2.0 notice is retained in
[LICENSE.kickmix.txt](LICENSE.kickmix.txt). See [REUSE.md](REUSE.md).
