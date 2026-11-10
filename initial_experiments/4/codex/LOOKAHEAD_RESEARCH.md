# Look-ahead moment maps and in-place series blocks

These experiments are excluded from the archive. The current selected score is
9.367722153134604. Experiments are stopped at the user's request; the 14.0
target remains unmet. All described jobs are terminal. See [EXPERIMENT_REPORT.md](EXPERIMENT_REPORT.md).
The accepted new DQI kernel is documented separately inCUBIC_TRACE_DQI_RESEARCH.md.

## A pivot-robust four-moment quantum component

For kept e,d and live a,b define D2=d^2+e*a and

    u=(d*b+a^2)/D2, v=(e*b+d*a)/D2.

The original a is recovered as e*u+d*v. `lookahead_pair_quantum.py` implements
this map in place when D2 is nonzero, including e=0. It selects e as pivot when
nonzero and d otherwise. For e=d=0 and singular D2 it supplies an explicit
reversible extension rather than applying an invalid zero-scaling erasure.
Both temporary zero flags are repaid before return.

The quantum test passes16384shots in each forward, inverse and one-builder
composition direction, including forced e=0, d=0, D2=0, and e=d=0 cases,
exact source/output values, clean scratch and zero phase. It uses102Q including
48public qubits,654Toffolis per direction. Evidence:
`results/lookahead_pair_quantum_check.json`.
This is a total reversible component with Padé meaning on its stated nonzero
minor domain; it is not yet a complete Euclid block or decoder gain.

## Six-moment closed-form coordinates

`lookahead_moment_map.py` maps F=(e,d,a,b,c,f) to the two adjacent monic Padé
polynomials B2=(u0,u1,1), B3=(v0,v1,v2,1), and reciprocal normalizer. Put

    D2=d^2+e*a,
    u0=(d*b+a^2)/D2, u1=(e*b+d*a)/D2,
    p=c+a*u0+b*u1,
    alpha=(f+d*u0^2+b*u1^2)/p, z=p/D2,
    v0=alpha*u0+z*d,
    v1=u0+alpha*u1+z*e, v2=u1+alpha, pinv=1/p.

This requires D2 and p nonzero but allows e=0, which breaks an ordinary
first degree-one quotient. Its inverse is explicit. Let
r0=v0+alpha*u0 and r1=v1+u0+alpha*u1. Then

    z=(r0^2+r1^2*u0+r0*r1*u1)*pinv,
    e=r1/z, d=r0/z,
    a=e*u0+d*u1, b=d*u0+a*u1,
    c=1/pinv+a*u0+b*u1,
    f=alpha/pinv+d*u0^2+b*u1^2.

For a larger cofactor block the transfer rows are

    A2=(d+e*u1,e), B2=(u0,u1,1),
    A3=(a+d*v2+e*v1,d+e*v2,e), B3=(v0,v1,v2,1).

Their determinant A2*B3+A3*B2 is the constant1/pinv. The12000-case classical
probe accepts8994nondegenerate cases and checks the inverse, full Euclid
cofactors, and this determinant. It includes3002paths with quotient degrees2,1,
and5992ordinary1,1,1paths. All3000forced first-two-zero cases are rejected by
D2=0, along with6other degeneracies. Evidence:
`results/lookahead_six_moment_map.json`.
The six-field quantum implementation below now exists. A full global cofactor
schedule using those coordinates has not been built.

## Reversible six-moment implementation and opening block

`lookahead_norm_reciprocal.py` implements the involution

    (e,d) -> p*(e,d)/(d²+e*d*v+e²*u),

retaining u,v and nonzero p. It chooses a nonzero coordinate before forming
the ratio t=d/e, computes H=t²+t*v+u, and sends e to p/(e*H). A zero H is
armed to one; the zero pair stays zero. Both temporary flags are cleared.
4096 cases per forward/inverse/composition direction pass, including forced
zero norm and zero pair, with clean scratch and phase. Q101 including60public,
854T per direction; `results/lookahead_norm_reciprocal_check.json`.

`lookahead_six_quantum.py` combines this with the prior two-moment solve. It
stores (v1,v0,u0,u1,pinv,v2) in the original six physical words and supplies a
total reversible extension at singular D2 or p. The nonsingular branch agrees
with the independent six-moment formulas, including e=0. Polynomial
identities reduce the p and alpha shears to four field products each; the p
inverse is shared with norm scaling. Q126 including72public,2094T per direction.
16384 cases pass each forward/inverse/composition direction, including forced
singular cases, zero inputs, scratch and phase.7015 cases also pass the
independent nonsingular Padé reference. Evidence:
`results/lookahead_six_quantum_check.json`.

`opening_euclid_lookahead.py` applies the map to the initial six syndrome
moments and updates the remainder by multiplying the series by reversed B3,
preserving its first six original coefficients until their conversion. A
temporary six-word bank permits reversible copy/transform/erase. N22 passes
1024 cases each forward, inverse, forward-then-inverse and fresh arbitrary
inverse-then-forward composition, including zero and singular heads, scratch
and phase. Q390/T6078 per direction. Evidence:
`results/opening_lookahead_check.json`.

This opening-block construction loses the resource comparison. At N70 the
four-degree prefix, including one following serial step and a persistent zero
normalizer guard, costs28114T/Q970; the current first-step-plus-three-step
prefix costs20360T/Q943. This is a resource comparison only; the new prefix
has not undergone complete decoder verification and is excluded from the
archive. `results/opening_lookahead_prefix_profile.json` records the counts.
The corrected copy/erase direction is covered by the final passing tests.
No selected decoder dependency was changed by these experiments.

## Direct terminal and full-series blocks

`terminal_euclid_block.py` exploits a final block with no remainder tail. It
transforms the original remainder words directly into local quotients, then
consumes them into the cofactors once. There is no duplicate head bank or pair
of speculative cofactor passes. AtN22/b4 it passes512shots each forward,
inverse and fresh arbitrary-input composition, including zero remainders,
flags, scratch and phase: Q295/T8840per direction.

`series_euclid_block.py` generalizes the representation change. Convert the
whole retained remainder to its local ratio, perform b Schur quotients, update
the cofactors once, then convert the surviving ratio back to a remainder.
The exact identity passes1000independent classical blocks. AtN22/k7/b2 its
quantum forward/inverse/arbitrary composition tests pass512shots atQ293/T7024.

These schemes save workspace but lose gate count. All78tested N70full-series
blocks at k19..34 are more expensive than the same selected serial steps.
Even the terminal blocks lose slightly: final width1 costs4253T versus4163,
and width8 costs47758T versus46156. The full-series blocks use902Q with all
35flags allocated, but replacing late steps does not reduce the earlier961Q
peak. Evidence: `results/series_block_native_schedule_comparison.json` and
`results/terminal_block_n70_profiles.json`. No full decoder was built from
these nonwinning profiles.
