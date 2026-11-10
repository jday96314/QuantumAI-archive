# Current alternative architectures (research only)

The authoritative archive is now **9.367722153134604**, SHA256
`9168da25a9ac35f661885edad3e208db215bf67fa6525cca0ccc5d8831318b68`.
The accepted inverse change is inZERO_HEAD_AFFINE_RESEARCH.md.
The separate accepted field-kernel change is inCUBIC_TRACE_DQI_RESEARCH.md.
None of the work below changes its selected circuits or native evidence.
Experiments are stopped at the user's request; the 14.0 target remains unmet.
All jobs described here are terminal. See [EXPERIMENT_REPORT.md](EXPERIMENT_REPORT.md).

## Direct Boolean phase compilation

`boolean_polynomial_phase.py` represents the existing polynomial, affine
restoration and phase shift as a74-input XOR/AND graph. Its795ANDs match the
independent integer model on32857 boundary/random inputs. Carry-save
arithmetic instead gives863ANDs, also model-correct.

`spooky_phase.py` and `spooky_linear_phase.py` implement measured eviction
with classical phase ghosts and eventual repair. Both pass their random-network
quantum checks. The latter streams wide linear operands instead of pinning all
parents. Neither scheduler produces useful resources on the actual polynomial:
small budgets exceed the20000-recomputation screen; at512cached nodes the
first needs2323ANDs/Q589, and the second1808ANDs/Q595. Large caches have783live
AND computations but over850qubits. These results reject those greedy
schedulers, not all possible Boolean compilers. See their JSON profiles.

ABC resynthesis of the complete795-AND front end produces1236ANDs, worse than
its input. `results/boolean_poly_abc.log` and
`parameters/boolean_poly_mapped.blif` retain that experiment.

The measured-eviction idea is standard background:
[spooky pebbling](https://arxiv.org/abs/2110.08973),
[improved spooky pebbling](https://arxiv.org/abs/2401.10579).
No external implementation was copied for these local compilers.

## Seed/residual inverse square root

The follow-up `goldschmidt_row_model.py` now uses signed radix-four row
rounding for both products, rather than whole-product rounding. It forms
R=R0-Booth(V,z), then S=R+square_correction, and returns
Y*2^(outbits-seedbits)+Booth(Y,S). The square uses a small rounded or floored
high part t of R; monomial variants separately round its signed quadratic
terms. Every recorded numerical screen covers all3133452 inputs.

Examples with actual row rounding and a floored monomial square pass:
f8/s10/r27/o27/square18 has maximum relative error4.43824741403e-8;
f8/s10/r28/o27/square18 has4.37799936304e-8;
f7/s10/r29/o28/square18 has5.74285208277e-8.
The two f7/o27 rounded-square models fail. `results/gold_row_checkpoint.json`
collects13 final numerical models and the quantum component evidence.
No seed/residual model has been built into an accepted complete circuit.

`split_square_correction.py` supplies a reversible low-workspace square
correction. The high part of R stays fixed while its signed monomials update
the low part plus a signed three-bit counter; integrating that counter
temporarily exposes S. The inverse first unintegrates, then clears the counter
against the unchanged high source. Controlled3*2^k constants share one adder
instead of two increments. It passes16384 cases each forward, inverse and
both fresh compositions, arbitrary incoming counter, scratch and phase.
For a20-bit R split at9, square shift11 and37 additional live qubits, it uses
71Q/425T per direction under a90-qubit limit. The21-bit/split10/shift10 version
also passes. An earlier deliberately oversized ambient56 case could not fit
90Q; the final ambient37 result is the relevant staged allocation.
The square remains expensive:425T is only this component, before seed-table
replay, both residual products, outer products, phasing and normalization.
These counts do not establish a complete resource win. The fixed-top/clean-
headroom change discovered during this work is a separate inverse candidate.

`goldschmidt_model.py` independently evaluates a short seed and its exact square,
a centered residual, and `y0*(1+r/2+3*r*r/8)`. Every reported model checks all
3133452inputs. Passing examples include f7/s10/r25/o27, maximum relative error
5.2798712e-8, and f8/s10/r25/o25 with one combined correction rounding,
4.58907521e-8. The former uses256lookup rows and seed/square/residual-center
widths10/20/17. No quantum implementation or resource improvement follows:
these models round whole products, whereas the existing Booth emitter rounds
individual signed rows. All detailed models and tables usegoldschmidt prefixes.

## A compressed cubic quantum circuit

`fit_cubic_table.cpp` fits a cubic using the actual radix-four row-rounding
model. Coordinates usefrac6 and a finer7-bit lookup:256rows including parity.
Widths27/18/11/5 fail109leaves (worst8.40215436e-8).
Widths27/19/11/5 pass at5.64446021034e-8;
28/19/12/5 pass at4.59077713622e-8. These are complete numerical model checks.

`cubic_phase_invsqrt.py` implements the27/19/11/5 case with a new intermediate
cleanup schedule. After a Horner source is consumed, undo its high carry
integration, measure its now lookup-only high bits, then restore those bits and
integration before the inverse arithmetic. It also borrows conditionally clean
normalized-input bits for phase headroom, using the previously checked helper.

The first compiled version exposed a23-bit affine-range overflow. The safe
version uses the finer7-bit address for affine preconditioning and a2^17bias:

    R0 = C0 + q_fine*2^18 - 11*2^23 - parity*2^25 - 2^17.

After the entire polynomial this accumulator ranges16318..8257535 over all
inputs, within23bits. The one leaf that rounded the normalized value to exactly
1 has its constant lowered by one unit so the27-bit phase mantissa stays below
1. The all-input maximum relative error remains5.64446021034e-8.
The safe table is `parameters/cubic6f7_w27_19_11_5_safe.txt`.

Corrected candidate `candidates/cubic_phase_corrected_q94.kmb` passes131072
quantum shots against the independent integer model, including4096consecutive
small inputs, random large inputs, exact phase, source restoration and all
scratch. Evidence: `results/cubic_phase_corrected_q94_quantum_check.json`.
It uses94qubits/2458staticToffolis, so it loses to selected90Q/1816staticT even
before rotation synthesis. It is excluded from the archive. Earlier resource
screen files and uncorrected candidates are historical rejected versions;
current source defaults rebuild the safe corrected candidate.

Reproduce:

```bash
source env.sh
"$EXP_PYTHON" cubic_phase_invsqrt.py
"$EXP_PYTHON" check_cubic_phase.py --shots 131072
```

## Direct Walsh phase synthesis

`walsh_invsqrt_probe.py` transforms all2^21normalized mantissa values and keeps
the largest Walsh coefficients. At1024terms the maximum relative error is
3.92696497e-7;2048terms reach3.71410814e-8. This favorable screen excludes parity
and exponent modulation and still needs2047nonconstant rotations. It establishes
numerical feasibility for this truncation, not a winning circuit. Evidence:
`results/walsh_mantissa_21_probe.json`. No native circuit was built.

## Fixed sparse-dictionary normalization certificate

`fixed_dictionary_bound.py` considers all76one-hot vectors and2850pair sums.
Their packed projector feature matrix P is invertible. If H is its coefficient
covariance, then for every factor-column gauge

    sum_c ||w_c||_1^2 >= sum_ij |H_ij| >= <S,H>

for any fixed entrywise sign matrix S. Let B=P^-T, d be the packed identity,
and A=B^T S B-d d^T. For any physical approximation within epsilon in Frobenius
norm, the native normalization algebra gives

    u >= -c + (<A,Gref>-epsilon*||A||F)/2,
    lambda_eff^2 > gap_floor*u.

Taking S=sign(B Gref B^T), epsilon20mHa gives a floating-point certificate
u>=992.113769 and lambda_eff>=62.321940. The feature inverse round trip has
maximum error1.11e-16. At the tested brot3/C276,C300,C350 shapes the heuristic
alone costs1.824b,1.976b,2.293b respectively, versus selected708.547m including
angles. These are specific resource screens; the normalization inequality is
independent of rank and column gauge for that fixed full dictionary. They do
not bound other dictionaries or the overall competition score. Evidence:
`results/full_pair_dictionary_bound.json`.
