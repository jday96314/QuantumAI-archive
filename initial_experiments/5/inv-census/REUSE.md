# Reuse and attribution

Read-only origin: `initial_experiments/5/inv-width/`.

Copied Python modules: binary_centered_invsqrt, affine_invsqrt,
fixed_top_affine_invsqrt, zero_head_affine_invsqrt, booth_adder, short_adder,
monomial_unlookup, optimize_workspace, native_sample. Their original docstrings
and algorithms are retained. `compact_wires.cpp` is the prior wire-lifetime
compactor. The prior coefficient table and accepted circuit/evidence were
copied without changes. These prior artifacts credit round 4 and earlier utility
sources; the host kickmix package is Apache-2.0, copyright its respective authors,
and its original notice is copied as LICENSE.kickmix.txt.

Local changes: binary builder adds headroom-free exponent phase modes and
allocates delayed-lookup scratch on demand; affine builder permits zero permanent
phase headroom and handles a zero deferred allocation without slicing away its
whole register. Experimental affine restore allocates private temporary carries.
No copied module or host dependency outside this directory was changed.

All other census kernels, exact certificates, models, checks and reports were
written here. No source from https://github.com/andigu/uncompute was copied or
imported. Gu et al.'s numerical architecture ideas are cited in PHASE1.md and
RESULTS.md; this is conceptual reuse, not a copied implementation.
