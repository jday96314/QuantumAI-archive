# Affine inverse-square-root coefficient storage

The current selection is now `candidates/zero_head_affine.kmb`, fully native
verified at182956.5271606425cost,90Q/1793staticT. The newer fixed-top and measured
headroom changes are documented inZERO_HEAD_AFFINE_RESEARCH.md. The historical
affine circuit and all its evidence below are retained unchanged.

The selected circuit is `candidates/affine_best.kmb`. The full unchanged host
verifier passes **200,540,928 shots**, with cost **185,026.55749858925**, 90
qubits, 1805.0004279625 mean Toffolis, 1816 static Toffolis, and 34 rotations.
The worst error/tolerance ratio remains 0.97967143955 at x=958102.
`results/affine_verify.json` records the result and native circuit hash.

Write q for the seven-bit coarse interval index and p for normalization parity.
The original interpolation center is `a=1+(q+1/2)/128`. Keep the existing slope,
curvature, and rounding rules. Replace the stored constant coefficient C0 by

    R0 = C0 + a/4 - A[p],  A = (961/1024, 1217/1024).

All quantities here use the original 27 fractional bits. The transformed table
constant fits 23 bits. Horner evaluation consequently gives the original
polynomial value plus the same interval-dependent affine offset. Before phasing,
restore `A[p]-a/4` modulo the 27-bit core.

The full classical all-x residual range is147391..8388606, inside0..8388607;
see `results/affine_center_only_ranges.json`. The upper bound is only two integer
units below overflow, so changing these offsets or the table requires checking
the range again.

The normalized input's top bit supplies 1/2 in the phase word. The remaining
restoration is `3/16 + p/4 - q/512`: two constant increments, a parity-controlled
increment, and one seven-bit address subtraction. Reverse these operations after
the phase. Three extra phase-headroom wires allocate only after high curvature
bits have been released. The selected schedule also loads residual bits21..22
late and measures them early, with joint lookup-phase repair.

`affine_invsqrt.py` applies checked source replacements to a copy of the existing
builder; it does not mutate the production quadratic builder or host code.
`parameters/affine_build_source.py` is the most recently generated copy, and may
reflect a research schedule. `build_affine_invsqrt.py` pins the selected schedule:
limit90, keep4, delay2, product17, deferred3. Its rebuilt circuit matched
SHA256 `89dc5547e068941b138fab5a1d68d33877ab43230a1b180043a84282852d34c1`
byte-for-byte before packaging.

The prior 9.2155702694 archive, summary, and audit are preserved under
`revisions/delayed/`. This inverse change reached **9.243872873372679** overall;
the later chemistry rotation reaches9.243907050802967. The14.0 goal remains unmet.

Other fitted coordinate widths were screened with actual compacted wire counts
in `screen_quadratic_fits.py`. Their smaller coefficient registers did not offset
larger lookup or arithmetic costs. Fine8 fitted tables pass the full classical
all-x model; these model checks are not native submission verification. An
earlier signed-slope affine prototype also passed a small smoke check but was
more expensive. Only the center-offset version above is selected.
