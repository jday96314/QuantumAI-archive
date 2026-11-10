# Compact decoder blocks and completion padding

These are research components, excluded from the archive. The current fully
verified score is **9.285435820940434** after the separate angle-gauge change.
The goal remains ACTIVE at14.0.
No complete replacement DQI circuit or score improvement follows from these
tests. All new files are in4/codex; upstream source remains read-only.

## Implemented arithmetic and block map

`block_euclid_ratio.py` checks2590 classical quotient blocks at multiple stages
of80random N70 instances. With `P B+R Bprev=p*x^N`, set T=reverse(R),
A=reverse(B), D=reverse(Bprev). The first2b coefficients of

    F = T*A / (p + z^2*T*D)

determine b normalized quotients `(alpha,beta)` independently of global degree.
Degree-gap blocks are counted and skipped by this diagnostic.

`mobius_series.py` implements this map without assuming T[0] nonzero. Conjugate
one unit-series reciprocal by Toeplitz maps and field scaling. Both forward,
inverse, one-builder composition, phase, and scratch pass512shots, including
zero T and zero T[0]. At n4,k3 it uses177Q,900T/direction; at n16,k8 it uses441Q,
9048T/direction. Both have33extra qubits beyond their research public operands.
`p_is_inverse=True` lets the compact decoder supply pinv directly and is covered
by the block tests. `BoundedArithmetic` avoids out-of-range coefficients and
chooses direct pairs when cheaper than Toom. Original unoptimized checks were
preserved with `_original_check.json` names.

`monic_euclid_step.py` implements a complete single generic compressed step,
assuming nonzero pivots. N20,k3 passes512shots at285Q/3354T per direction;
N70,k17 passes at885Q/15030T. Both check independent classical outputs, inverse,
one-builder composition, phase and scratch. `monic_euclid_local.py` composes
the first step and generic steps into a local quotient transform: N8 uses
141Q/2376T, N16 uses237Q/12306T. Each passes512shots with gaps explicitly skipped.

`continued_fraction_series.py` instead stores each quotient pair directly:
beta=F[0], alpha=(beta/F)[1], next F=(beta/F)[2:]. Unit reciprocal never reads
the physical F[0] word, which retains beta. Optional retained flags arm zero
beta to1, with record/replay modes for duplicate transforms. Arming preserves
reversibility but alone does not preserve a completed decoder's locator.

`blocked_monic_euclid.py` copies only the2b head words to a temporary bank,
extracts quotient pairs, and makes three cofactor passes. With Bfinal known,
apply its unit upper Toeplitz map to the full C array and undo that map only
on the head. Restore old cofactors to divide the remaining tail by old B.
Finally transform the restored head to the same quotient pairs, XOR it to
zero, consume the copied pairs into the final cofactors, and swap those into
the now-zero public head. This avoids a copy of the old full cofactor bank.

The exact tail identity on a generic block is

    reverse(Cnew) = ((reverse(Bfinal)*reverse(Cold))[2b:]) / reverse(Bold).

Discard the top2b numerator coefficients **before** the old-cofactor division.
Commuting the division past this truncation gives an incorrect formula.

| Quantum block | Q | T/direction | Shots |
| --- | ---: | ---: | ---: |
| N20,k3,b2 |321|7692|512 PASS |
| N20,k3,b2, zero arming |323|7714|512 PASS |
| N70,k17,b4 |969|40590|512 PASS |
| N70,k17,b8 |1065|95298|512 PASS |

Each includes classical values, inverse, composition, scratch and phase.
The armed test includes all-zero remainders and zero leading coefficients.

`blocked_euclid_sequence.py` composes the first armed step and all blocks,
retaining one flag per quotient. N20,b4 passes512shots, including all weights,
zero syndromes and forced gap inputs, at379Q/35612T per direction. Its output is
a reversible compact state, **not a completed decoder**. Generic fake padding
can corrupt the locator after early completion; the next construction fixes it.

## Completion, locator and magnitude

`check_compact_forney.py` checks the unpadded compact state at every nonzero
weight1..35:698cases,12548error roots PASS;2degree gaps skipped. B is the monic
locator. For H^2=(B Bprev)' and an error root x,

    error_value = (1/pinv)*x^69 / H(x)^2.

Scaling H by sqrt(pinv) lets the50T inverse-square phase kernel use only H(x)
plus the known classical position power. This is a classical identity check.

`coprime_padding.py` keeps completed paths at full degree without introducing
field roots. On completed quotient pairs, use distinct irreducible quadratics
P=x^2+x+b (absoluteTrace(b)=1):

    Bnew=P*B, Bprev_new=P*Bprev+B, Hnew=P*H+B.

At most one unpaired completed quotient uses P=x with omega=3362 in GF4:

    Bnew=x*B
    Bprev_new=omega*x*Bprev + omega^2*B
    Hnew=omega^2*x*H + omega*B.

Pairs start at global stage1, so an unpaired completion can occur only once.
All artificial quadratic factors are distinct, and x has no root among the
nonzero error positions. Thus no added factor overlaps a true locator root.
After a fixed-degree Wronskian map, synthetic division unpads B and H using
only affine field operations; apply sqrt(pinv) scaling **after** this unpadding.
This keeps the unpadding linear in B,H.

The original all-linear-factor probe passed699cases but had28artificial-root
collisions. The quadratic variant's larger check passes **34850cases and
626265error roots**, all weights1..35, with **zero artificial collisions**;
150generic degree gaps are explicitly skipped. Evidence:
`results/coprime_padding_n70_b4_quadratic_check.json` and `_large_check.log`.

`quantum_coprime_padding.py` implements controlled affine extensions of both
maps. Extra quotient slots carry polynomial remainders, making each map a
bijection on arbitrary input words. At k17, quadratic padding uses445public Q,
1569T, and linear padding421public Q,2385T. There is **no extra arithmetic
workspace**. Both pass512shots for forward, inverse, composition and phase,
including random remainder words and both control values.

The `padding=True` option in `blocked_monic_euclid.apply` integrates the maps
in pairs. First do the ordinary cofactor updates. A controlled linear inverse
of the beta=1,alpha=0 specialization restores the old pair on a completed
branch; then controlled padding replaces it. A joint flag chooses quadratic
padding, and second-flag AND NOT first chooses the sole linear pad. All other
inputs still have a literal inverse, including genuine degree gaps.

`check_padded_sequence.py --n 22 --block 4` passes512
shots on every supported weight, plus zero syndromes, against independent
polynomial algebra. Its composition test uses fresh arbitrary words, forced
zero leaders and zero syndromes. Flags, scratch and phase restore exactly.
Resources: **406Q,65605T per direction**. Evidence:
`results/padded_sequence_n22_b4_check.json`. Wronskian and decoder phase remain
unimplemented in this sequence.

## Resource evidence and next architectural work

`profile_blocked_euclid.py` counts normalized emitter product terms without
materializing the huge Clifford tape. Direct controlled padding gates are
added from the builder's actual gate counter. Its N22 padded whole-transform
count **exactly matches406Q/65605T** from the compiled quantum check.

The generic mixed schedule in `results/blocked_euclid_schedules_n70.json` gets
278370T/direction at969Q before first step, zero flags, completion, H, scan and
repetition. It is an arithmetic-only schedule, not a compiled full decoder.
It mixes serial late steps with early blocks, and cannot simply be substituted
into the current even-pair completion implementation.

Whole forward N70 resource counts (including first step and all35zero flags):

| Block | Unpadded Q/T | Completion-padded Q/T |
| --- | --- | --- |
|2|956 /391609|958 /667729|
|4|1004 /332017|1006 /608137|
|6|1052 /343291|1054 /619411|
|8|1100 /369505|1102 /645625|

These exclude Wronskian, unpadding, scan and repetition. **The current padded
implementation is too expensive to establish a DQI score win**. Padding adds
276120T per direction across all block sizes because all three cofactor passes
pay its controlled linear maps. Do not build a huge whole submission on the
assumption that the low generic estimate already includes completion.

Potential follow-up: apply cheaper monomial padding during the decoder, unshift
the pair after it, and add coprime factors only once before a fixed-degree H
map. Or fuse/block the cofactor matrix update rather than paying O(k*b) field
scalings in each of three passes. These are proposals, not tested savings.
The true weighted distribution has P(L<24)=.00184553 and P(L<25)=.00507479;
late-only expensive handling is another measurable tradeoff, not a correctness
claim. Native verification must price every accepted phase failure.

For the fixed-degree H map, reverse both monic polynomials C=rev35(B),
D=rev34(Bprev). Split C=E^2+z O^2, D=U^2+z V^2. Then

    reverse34(H)=E*U+z*O*V.

E[0]=U[0]=1. This matches the existing generic Sylvester prototype with a
constant payload word, avoiding a new map derivation. Two temporary constant
words fit below the decoder's peak. A constant translation x->x+1 of B,Bprev
before this map and H afterward can remove the deliberate zero root from the
source's highest-degree slot; translation is linear and Clifford-only.
**This adapter has not been built or tested.** The current generic Sylvester
prototype also needs zero-pivot/value-restoration audits before integration;
its prior tests skipped degree gaps. Inspect psi_forward's existing epsilon
repair before changing it: the original generator does handle degree defects,
while the payload scaling/quotient modifications impose additional assumptions.

Read-only arithmetic background consulted: Dumas and Grenet,
[Fast in-place accumulation](https://membres-ljk.imag.fr/Bruno.Grenet/publis/DuGre26.pdf).
The implemented block identities, circuit schedule and padding probes above
are local experiments; the paper is not evidence for their quantum correctness.
