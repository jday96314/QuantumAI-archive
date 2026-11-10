# Read-only sources and attribution

The unchanged host DFTHC verifier/math/circuit infrastructure is imported from
`../../../from_host/verifiers/dfthcsos` (Apache-2.0; license in
`../../../from_host/verifiers/LICENSE.txt`). No host code or reference inputs
are modified.

`core.py`, `psd_capacity.py` and `verify_certificates.py` adapt QR, covariance,
rank-refit and dual-certificate methods from completed `../chem-psd` and
`../chem-discrete`. `descend.py` adapts the completed `../chem-qfit/discrete.py`
physical gradient/projected-curvature proposal method for arbitrary trees.
`physical.py` adapts completed `../chem-qfit/physical.py` offset, one-body and
column-gauge calculations. These source experiments are read-only.

`oracles.py` adapts completed `../chem-qfit/oracles.py` and imports read-only
round-4 Codex `build_selectswap_chem`, `monomial_unlookup` and `jacobi_l1`.
The shared low grid-phase-bit extension is local to this experiment and retains
the host's integer output semantics and accounting.

All new tree generation, encoding, binary-partition equivalence, topology
search/replay, local binary restructuring and reporting code is local.
`plot_results.py` uses the completed `../chem-psd/.plotdeps` matplotlib
installation read-only; caches are directed here. No dependencies were installed
or shared environments changed.
