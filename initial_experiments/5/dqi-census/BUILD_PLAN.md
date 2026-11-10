# Conditional build plans and integration gates

There is no component-supported whole decoder at >=2x. Phase 3's measured
projection gate is not met. The following is a concrete staged research/build
plan for the still-open resource targets; it must not be described as a
ready-to-build complete decoder. The unresolved primitive is identified in
step 2 rather than hidden in a zero-cost adapter.

## First priority: resident rational coordinates, C

Resource target: Q940; core240000T including forward and reverse; adapters10000T;
80T*4024.017333984375 scan =321921.38671875T. Total571921.38671875T.
Using the unmodified host CP at zero failures, cost602003162953.9934 (2.6548x).
At15 failures cost627940177452.299 (2.545x); at24 cost638647188678.255, barely
inside the2.5x target639281691753.9706. Use <=15 failures per full validation
stream as the engineering margin, not as a changed verifier rule.

1. **Specify the resident prefix permutation before gates.** Work in two-symbol
   steps. On a generic length2k prefix, normalized sigma has k free coefficients,
   omega has k coefficients, and together they recover the prefix by
   S=omega/sigma mod x^(2k). Logical slots are:

   | Physical slots | Logical contents |
   |---|---|
   | S[0:k] | omega[0:k] |
   | S[k:2k] | sigma[1:k+1], sigma0 implicit1 |
   | S[2k:70] | untouched future syndromes |

   At an even boundary this already fills ALL processed slots. A full tau
   template does not fit there as additional data. It must be recovered a
   coefficient at a time from the rational state, or encoded by a proved small
   descriptor. Specify D_(2k) on zero pivots, rank defects, zero syndromes and
   completion; retain remaining degrees of freedom in the unused padded
   coefficient positions. A rank tag alone cannot erase data. Define the pair
   update explicitly as D_(2k+2) (D_(2k)^-1 x I_2). This defines its required
   permutation, not a cheap implementation. Check its inverse independently.

2. **Replace the failed disjoint update by a real resident primitive.** Implement
   at k=1,2,8,17,34. It must accept the slots above, two new syndrome words,
   and <=100 clean temporary qubits; materialize no separate length-k sigma,
   omega or tau bank. At k34, the last two fresh syndromes finish the840Q
   allocation rather than providing extra clean memory. Stream a correction
   descriptor with an explicit recovery formula and cleanup proof. The current
   reference V -> (sigma_new/sigma_old)V computes correctly but requires a
   disjoint sigma/tau and costs31962T at k35. Its70-update realization costs
   at least4.8066mT before correction generation. Sharing its sigma with a
   zero V row is invalid: that row holds sigma, not logical zero, when the
   triangular emitter reads it. A new alias-safe identity is required.

   Gate target averaged over35 pair updates, both directions: <=6857T/pair.
   The current degree35 reference costs63924T for one forward/reverse pair,
   ~9.3 times that average allowance. Coefficient-zero trimming alone is not
   an implementation of the required algorithm.

3. **Prove the branch/history representation.** Reserve at most20 persistent
   history bits at the scan and an80Q scan scratch allowance (840+20+80=940).
   A scalar discrepancy may need12 bits; degree/defect tags need the balance.
   If more are required, reprice Q immediately. Do not discard a template on
   the grounds that it is a deterministic function: show and charge its
   recovery circuit. Reuse dqi-dispatch RESULTS as algorithmic evidence only;
   do not import its changing builders. A full state comparison against two
   independent classical decoders must cover arbitrary-stage gaps, consecutive
   gaps, genuine early completion, and the final boundary.

4. **Emit a small complete oracle, then N70.** At N8/N16 test native zero, one,
   maximum-weight and forced degree-gap cases with independent per-shot masks;
   check syndromes, every scratch bit, and exact phase. The phase midpoint
   must hold sigma and omega in the declared physical slots. Insert the
   already-validated80T point kernel; port the70-point classical interpolation
   and reverse all mask maps. Explicitly normalize omega consistently with
   sigma. Count all input/control paths as emitted, irrespective of rarity.
   If width exceeds940 or core+adapter exceeds285415T at940Q, the2.5x
   zero-failure budget fails; rescore without pretending a new bound rules
   out the architecture.

5. **Final acceptance only after complete native integration.** Preserve the
   incumbent until a complete N70/M4095 circuit passes default102400 shots
   and an independent salt. Commands (from this directory):

   ```bash
   source env.sh
   flock /tmp/r6-fullverify.lock "$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py candidates/resident_complete.kmb --max-workers 2 > results/resident_default.json
   flock /tmp/r6-fullverify.lock "$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py candidates/resident_complete.kmb --max-workers 2 --salt dqi-census-independent-a > results/resident_independent_a.json
   ```

   These filenames are future outputs; neither circuit nor full-run reports
   exist. Salt is predeclared, not searched. Require exact source/artifact
   hashes, unchanged host confidence, and the actual R Q² mean-T in both
   reports. Promotion/packaging is separate from this component census.

## Second priority: direct F' converter, A

A useful measured base case exists. For B=x²+b1*x+b0 and A=a0+a1*x:

```
u=a1*b1
N=a0²+a0*u+a1²*b0
D=b1*N
F'=(a0+u+a1*x)/D
```

The component computes F' by one scalar inverse in each XOR evaluator, then
uses XOR/swap/XOR to overwrite A while retaining an old-A work bank. This
is an involution on the pair even for singular data, so the full paired
inverse restores everything. On promised split squarefree B and unit A,
the work bank is zero at the phase midpoint.

The low-width emitter uses a parity specification for a1², releases N before
scalar inversion, and recomputes N for cleanup. One overwrite:960T,117Q;
existing A/B contribute48Q, so extra69Q comprises24 retained work+45 temporary.
Forward and arbitrary-state roundtrip each pass8192 native shots. B1 takes
all4096 field values in roundtrip controls, including the repeated-root case.
The faster792T variant needs81 extraQ, so it misses the74Q converter budget.
Neither a scalar-inverse shortcut nor this quadratic norm supplies a recursive
norm on an arbitrary degree35 split polynomial. Such a polynomial has no
known fixed quadratic tower; factoring/ordering its quantum roots is charged
work. Degree2 counts cannot be multiplied by17.5 and called a degree35 estimate.

Next steps:

1. Derive a total in-place modular inverse or direct reciprocal-product map
   on degree4/8 with <=74 extraQ, using the quadratic component as a base
   only where a proven reduction exists. An EEA/dialog alternative must
   reserve its actual live operands and payload. Khattar's4n² leading field
   products for direct division give205800T at n35 with42T products, before
   payload cleanup; that general scheme is not automatically within157k.
2. Profile BOTH directions at degree35, including old-A erasure, derivative
   preparation, padding/degree routing and singular branches. The fast
   Frobenius implementation already fails:1799952T per direction and11751Q.
   No additional full decoder build is warranted for that implementation.
3. Fold beta^70 into dynamic trace masks; insert11T scans and reverse the
   converter. Independent polynomial checks already cover288 cases across
   all degrees0..35,5040 roots and1179360 independent mask words.
4. At unchanged Q961/core649802, even FREE conversion only gives2.0639x;
   a real2x converter allowance is11248T per direction, much stricter than
   the157k threshold for merely beating the existing scan. With trimmed
   core624206 and projected Q948 counter layout, a30000T/direction map gives
  738262T and790376520895.1322 cost (2.022x), but new gap correction G is
   still unpriced. Exact target from budget.json is authoritative. This
   target has little slack and is not a2.5x plan.
5. If the complete projection remains>=2x including G and measured Q, use
   the same complete-control/default+salt acceptance flow above. Otherwise
   retain the converter as a component; do not select it by rare-branch
   frequency or a favorable salt.
