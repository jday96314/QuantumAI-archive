# Current selected phase kernel

The compact core below is unchanged. The audited score is now9.367722153134604
with the new88TGF16cubic-trace inverse-product scan. Full native102400shots,
independent full rebuild and both archive audits PASS. See
CUBIC_TRACE_DQI_RESEARCH.md for exact resources, hashes, math and rebuild commands.
The older checkpoints in this file are historical.

# Deferred squarefree padding for compact Euclid

The complete compact inverse-product decoder is now selected in the submission.
Current fully verified archive score is **9.310120488836658**; selected chemistry
remains `candidates/norm_cycle2_signed`. Experiments are stopped at the user's
request; the 14.0 target remains unmet. See [EXPERIMENT_REPORT.md](EXPERIMENT_REPORT.md). Earlier
completion and Wronskian prototypes below are research history.

## Selected complete decoder

`candidates/compact_product_slim_n70_m4095.kmb` passes the full native102400-shot
verifier with zero value errors and365phase errors. Q961, T1036334, native
repetition1.6820375683380908, final DQI cost1609838144370.2598. This barely
improves the inherited1612737419463.0198. Raw width savings do not imply a large
score gain: the generic-path degree gaps are now the material limitation.
The first assembled version passed but cost1678643393596.0247, so was rejected.

The selected version combines these changes:

* Monomial completion through the compact core; block3 through stage19, then
  single steps from stage22, with preserved-prefix tail multiplication.
* `controlled_linear.py` combines each fake-update correction and monomial
  permutation into controlled transvections, each oneToffoli after a Clifford
  basis change. Random affine maps, inverses and composition pass512shots.
* `low_workspace_inverse.py` releases norm-chain temporaries early using
  `Norm(x^-1)=Norm(x)^-1`, then measures inverse output before allocating its
  cleanup chain. Exhaustive4096-field-value/random-target test passes at the
  same102T. The helper explicitly requires an unjournaled emitter.
* Completion flags allocate as needed and release during the inverse.
* The scale word itself supplies the prior cofactor's leading coefficient.
  The locator's leading coefficient stays virtual: unary completion flags
  supply its value and derivative during evaluation by CNOTs.
* The inverse-product scan uses92T/point and avoids a Wronskian map and final
  squarefree padding. Both cofactor banks are rotated down by their known
  common monomial degree before evaluation.

The selected complete builder was executed again independently and its KMB
hash matches exactly; `results/compact_product_slim_rebuild_check.json` records
SHA256781e1c04fe1bca206adbf9f105ef537b31329dfbe53225c6bcbc4a07e3c2b2b4.
Both final ZIP audits pass. The preceding archive remains inrevisions/norm_cycle2/.
Rebuild commands are inREADME.md and the selected source manifest is
`parameters/compact_product_slim_source_manifest.json`.

## Direct replacement and preserved-prefix tail multiplication

The next two independently tested changes improve the same completion core:

* `direct_repad` replaces common x pairs directly with distinct quadratics by
  the upper triangular map `1+x^-1+b*x^-2`. It removes the separate cofactor
  unrolling and repadding. N22/block4 passes512 native shots in each direction
  and arbitrary-input composition, Q404/T49210 per direction.
* `BoundedArithmetic.triangular_suffix` multiplies the retained tail while
  preserving the copied head. This avoids computing and then undoing the
  head's triangular multiplication. With direct repadding, N22/block4 passes
  the same512-shot independent value, flags, scratch and phase checks at
  Q404/T46228 per direction. See
  `results/deferred_padding_direct_suffix_n22_b4_quantum.json`.

Full N70 forward exact counts, excluding H construction, H division, scan,
and repetition, are now:

| Block | Q | Direct repad T | Direct plus suffix T |
|---|---:|---:|---:|
|2|956|483326|476900|
|4|1004|423734|413318|
|6|1052|435008|421316|
|8|1100|461222|446480|

These are component counts and do not establish a new DQI score. The original
451699 block4 count below is historical. Reproduce the newer version with
`check_deferred_padding.py --n 22 --block 4 --trials 5 --shots 512 --direct --suffix`.
The exact counter is `profile_blocked_euclid.whole(70,4,deferred=True,direct=True,suffix=True)`.

A separate **rejected** `thin=True` experiment shortens the two temporary
cofactor passes to leading bands. It passed the valid-error polynomial values
but failed arbitrary-input composition: an independent512-shot rerun found
shot199 changed both values and phase. The failure is retained in
`results/deferred_padding_thin_rejected.json`. Its smaller resource counts
must not be used as an accepted component. The default remains `thin=False`.

## Mixed single-step schedule and an alternative phase scan

`monomial_serial_euclid.py` adds a single-step variant for the final rounds.
Its zero-pivot arming uses inverse(pinv), shifts the sublead, and masks fake
input contributions to the retained remainder's boundary. It then applies the
same monomial completion correction. Valid decoding values and arbitrary-input
composition pass512 shots at N14; a mixed N22 schedule passes1024 shots. The
N70 block4-through-k17 then serial-from-k21 schedule passes512 shots in forward,
inverse and arbitrary-input composition at Q1004/T345260 per direction,
including final direct squarefree repadding. See
`results/deferred_padding_direct_suffix_serial21_n70_b4_quantum.json`.

`adaptive_euclid_profile.py` records exact stage counts and width-constrained
dynamic programming schedules. At Q980, block3 through k19, then serial from
k22, costs362780 per direction including repadding. At Q1004 the mixed schedule
costs345260. These are still core counts, not complete decoder scores.

An alternative complete scan now avoids both squarefree repadding and the
Wronskian map. `compact_monomial_unpad.py` exposes the two variable-degree
cofactors by rotating their explicit monic coefficient banks down by the
completion degree. This costs4812T, with a linear XOR unary-to-binary map for
the flags; all degrees and the inverse pass512 shots. The prior cofactor is
scaled by pinv before evaluation. For native syndromes starting at index1,
the error value at root x is `x^N/(pinv*Bprev(x)*B'(x))`.
`inverse_product_kernel.py` supplies the92T root-gated inverse-product parity
kernel;4096 field cases pass sources, scratch and phase checks.

`compact_product_decoder.py` assembles this architecture. Its smaller N14,
m63 native verifier check passes2048 shots with zero value/phase errors;
this is not a competition result. Full-size build/verification work is tracked
in CONTINUATION.md. Do not infer a full-size score from the smaller test or
reuse the baseline repetition factor without measuring the new phase errors.

The previous completion schedule changed the cofactor pair during each of three
passes through every quotient block. Its full N70, block4 forward transform
cost608137T at1006Q, before the Wronskian, scan and repetition factor.

The new schedule in `deferred_euclid_padding.py` and the `monomial=True` option
of `blocked_monic_euclid.apply` separates two requirements. During the compact
core, common multiplication by x preserves the actual locator as a factor, so
the retained remainder stays zero after completion. At the end, divide both
cofactors by those known common powers of x, then multiply both by the same
squarefree padding polynomial. It consists of distinct irreducible quadratics
`x²+x+b`, with at most one x. The padded locator reaches degree35 and is
squarefree; its roots in GF(4096)* are precisely the actual error roots.

The cofactor pair may share the padding polynomial. The fixed-degree Sylvester
map needs the locator's even/odd halves to be coprime, which follows from the
locator being squarefree. It does not require this pair of cofactors to be
coprime. If their original Wronskian is H², multiplying both by P changes it
to(PH)². The known P can be removed from H afterward by controlled linear
polynomial division. That division and the H adapter are still unbuilt.

The first attempted deferred schedule simply undid the armed core's Fibonacci
dummy updates at the end. It failed when completion occurred inside a block:
later blocks could revive the retained remainder. It supported only499 of1100
classical cases. This was rejected, not counted as a successful generic handler;
see `results/deferred_fibonacci_undo_n22_b4_rejected.json`. Its quantum build was
stopped to correct the algorithm. The corrected core uses common monomial
padding in every pass, with the more expensive squarefree factors deferred.

Corrected evidence:

* N22, block4 classical:1098 cases,6585 roots,2 skipped generic degree gaps.
* N70, block4 classical:3486 cases,62665 roots,14 skipped generic degree gaps.
  Checks include exact dummy-factor removal, locator roots, squarefreeness,
  H multiplication/division, and inverse-square error magnitudes, all weights.
* N22 native quantum circuits:512 shots each for forward, inverse and composed
  directions. Forward/inverse values agree with independent polynomial algebra
  on every nonzero weight. The composition also uses fresh arbitrary syndrome
  words, zero syndromes and forced missing top pivots. Inputs, flags, scratch and
  phase all restore. Q404, T51675 per direction, versus the earlier padded
  prototype's406Q/65605T. This is a component check, not the host DQI verifier.
* The exact resource counter reproduces the compiled N22 count, recorded in
  `results/deferred_padding_profile_crosscheck.json`.

Full N70 forward resource counts from
`results/deferred_padding_whole_profiles_n70.json`:

| Block | Q | T per direction |
|---|---:|---:|
|2|956|511291|
|4|1004|451699|
|6|1052|462973|
|8|1100|489187|

These exclude H construction, division of H by the padding polynomial, the
root/magnitude phase scan and repetition pricing. The block4 core improves by
156438T/direction,25.72%, but its complete-decoder cost is unresolved. Do not
claim a new DQI score or build a huge full submission from this core count alone.
No research processes remain from these tests.

Reproduce:

```bash
source env.sh
"$EXP_PYTHON" check_deferred_padding.py --n 22 --block 4 --trials 100 --shots 512
"$EXP_PYTHON" check_deferred_padding.py --n 70 --block 4 --trials 100 --classical-only
"$EXP_PYTHON" profile_blocked_euclid.py --n 70 --deferred
```

Further arithmetic changes considered but not implemented:

* Three cofactor passes remain a major cost. Keeping an old cofactor long
  enough for the tail map, then erasing it from a short leading-coefficient
  reconstruction, could replace two full passes. Its memory and erasure cost
  need an actual derivation and profile.
* Nonmonic cofactors use unimodular two-term quotient multiply-adds, potentially
  avoiding per-word scaling. Their nonunit leading coefficients add scaling to
  tail transforms. Lazy normalization once per block requires reversibly
  clearing its scalar history; no complete cheaper schedule is derived.
* Remember that upstream `mon_slot` costs one field product(42T), with free
  measurement phase repair, not two products. `inv_word` adds102T for its
  complete lifetime, shared across every word scaled in its context. Earlier
  informal estimates treating each scaling as84T were too pessimistic.
* QD/Rutishauser and progressive QD recurrences were researched, including
  [NIST DLMF3.10](https://dlmf.nist.gov/3.10) and primary mathematical sources.
  Their rhombus rules alone do not establish an N-field reversible transform.
  An explicit invertible low-workspace schedule is still needed; no quantum
  QD circuit or useful resource count has been built.
