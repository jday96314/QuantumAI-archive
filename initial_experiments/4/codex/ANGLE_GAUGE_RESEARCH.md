# Exact signed angle gauges

The selected chemistry family is `candidates/signed_gauge7_prefix8`.
Full unchanged native verification passes **50880shots**, cost
**712929326.0737408**, lambda_eff48.89470055602427, total478Q/30504T.
The11individual oracles each use123Q/482T. All also pass2048random-target
shots each, restoring addresses, arbitrary incoming public words, scratch and
phase. Evidence: `results/signed_gauge7_prefix8_verify.json` and
`results/signed_gauge7_prefix8_random_targets.json`.

Backward Givens elimination normally uses positive hypot for each parent.
Instead retain the parent's sign in hypot and choose each rotation in
[-pi/2,pi/2). At a boundary, subtract pi and flip the retained parent sign.
The resulting preparation produces either u or -u. Both give the same
projector u*u^T and therefore the same one- and two-body operators.
`angle_gauge_probe.py` checks every vector up to global sign, maximum error
7.22e-16. The new8-bit words have bit7=bit6 **without exceptions**.

XOR differences of sign-extended words retain this equality. Before the
lookup, CX(public bit6, public bit7); load only bits0..6; then repeat the CX.
This writes the selected bit6 into both high bits while preserving arbitrary
incoming target states. It does not assume the target words are initially zero.

The select-copy difference bank shrinks from56 to49qubits. This permits an
eight-bit unary subtree for the nine-bit high address, rather than six-bit
subtrees. The public seven8-bit words retain56qubits. Peak ANGLE width falls
128->123, although the global chemistry temporary maximum is elsewhere, so
the scored478qubits remain unchanged. Total ANGLES cost drops11016.9725432
to10604T; the physical metrics and19900T heuristic contribution are unchanged.

`build_signed_gauge_chem.py --out candidates/signed_gauge7_prefix8 --lowbits 8`
rebuilds the family from the verified Jacobi candidate. The independent rebuild
`candidates/signed_gauge7_prefix8_rebuilt` matches all11circuit files byte-for-
byte and all7NPZ arrays exactly. Evidence:
`results/signed_gauge7_prefix8_rebuild_check.json`. The earlier full-unary
version `candidates/signed_gauge7` also full-verifies at714986032.7579296;
it is superseded by the prefix8 family.

The current archive's overall score is **9.285435820940434**. SHA256:
`5de660862765b7b4a2cefab452e1a731180f6c7a49e105fa30440feaca72e7ab`.
`results/archive_audit.json` binds all14members; the separate
`results/archive_independent_audit.json` reopens the ZIP, checks each member
against its native evidence and recomputes the chemistry canonical digest
`9db8ca09bffcd4ea9293e4b2eb3ef0b9c1962f26050300e4fd625d4e8f276883`.
The preceding9.2439070508 archive/summary/audit are in`revisions/jacobi/`.
The14.0 goal remains ACTIVE and unmet.

## Limits and other screens

The sign gauge alone does not make eight-word batches fit this select-copy
architecture: their795-row difference tables have binary rank56, so public64
+ difference56 + address10 already needs130qubits, before scratch. This is
a bound on this separate difference-bank layout, not on all possible oracles.
Unsigned canonicalization reduces early ranks but not late signed batches.

New inverse coefficient probes are outside the archive. `fit_adaptive.cpp`
has an optional `deep` mode that fully rescans the remaining candidates on a
failed smallest interval. Existing default behavior is unchanged; the separate
compiled binary is`fit_adaptive_deep`. Curvature precision6 at slope18 misses
5.9e-8, even with this local rescan; this is not a general impossibility proof.
Precision7/slope18 fits all-x but uses88Q/1886T at the best screen, worse than
the selected90Q/1816staticT circuit. Its fitted table has68different half-pairs,
versus13for the selected table, offsetting the smaller curvature register.

Slope19 with curvature8 or9 fits every coarse interval, and the builder drops
the unused fine address bit. The actual resource screens still lose: best
quick estimates are189980 and194028 respectively, versus185026selected cost.
Slope18/curvature9 is also worse at187740estimated. The widened inner slope
accumulator split is10, rather than9, to retain the original carry headroom.
All resource artifacts and exact all-x residual ranges are saved under
`results/precision_*_resources.json`; none was native-verified or packaged.
`affine_fitted_invsqrt.py` uses per-tag compaction temporaries. The initial
parallel screen had shared temporary names and was rerun after correcting
that race; only the corrected results are authoritative.

Read-only primary lookup background consulted:
[Motlagh and Pocrnic, Halving the cost of QROM](https://arxiv.org/html/2605.20334v1),
[Khattar and Gidney, conditionally clean ancillae](https://arxiv.org/abs/2407.17966).
The selected circuit's correctness and resources are established by its own
native and arbitrary-target checks above, not by those papers.
