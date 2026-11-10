# Small-ensemble chemistry factorization experiments

These classical fits are excluded from the submission. The verified archive
score remains9.304538097994794. The20mHa physical error limit is mandatory;
none of these alternative shapes reached it.

`grouped_thc_seed.py` agglomerates the35 reference ensembles by physical
covariance overlap. Each merged group retains its leading physical columns,
then selects a shared rank-one dictionary by simultaneous OMP.
`grouped_thc_fit.py` jointly optimizes normalized vectors and whitened factor
columns in float64 against the full reference matrix, including its negative
spectrum. The contracted loss is checked against direct matrix subtraction at
initialization. Trace normalization is optional; these runs enabled it.

| Shape R,B,C | Initial error | Adam10000 final best | Further L-BFGS1000 best |
|---|---:|---:|---:|
|10,64,35|3.078527Ha|0.130352166Ha|not run|
|16,64,24|2.492191Ha|0.063839810Ha|0.062566419Ha|
|16,64,64|0.059705702Ha after full-rank projection|not run|0.055777958Ha|

The full-rank projection holds the dictionary fixed and optimizes positive
semidefinite ensemble covariances without the original C24 restriction. It
shows a modest improvement for this initialization, rather than an acceptable
factorization. The L-BFGS runs lower the ridge from1e-8 to1e-10; their growing
weight norms and small physical improvement do not justify proceeding to angle
quantization or native circuit construction. All four training jobs are terminal.
These are results for the listed initializations, not global capacity bounds.

A separate low-precision sparse-coordinate dictionary probe uses one-hot
vectors and pair sums. Even selecting2500 of2926 pair-space coordinates leaves
38.812mHa physical error at rank350. At719 coordinates the error is505.905mHa.
The tested simple sparse dictionary does not provide a useful low-brot shortcut;
see `results/sparse_pair_projection.json`. Other dictionaries remain possible.

Resource-only shape estimates are in `results/grouped_shape_heuristic_costs.json`.
They exclude ANGLES and are not accepted native scores. Numerical artifacts and
logs retain the grouped10/grouped16 prefixes under parameters, candidates, and
results. Reproduction uses the torch2.13 Python with env.sh sourced.

## Later R24/R28 runs (all terminal)

The current verified archive is9.310120488836658, superseding the number above.
R24/B64/C32 initializes at1.963196460Ha. Adam12000 reaches24.645975575mHa;
L-BFGS1000 continues to22.206614165mHa, still invalid. Its weight squared norm
rises to625.47. Full-rank fixed-dictionary projection of an earlier snapshot
reaches24.642950521mHa. Alternating PSD covariances with dictionary Adam4000
reaches24.433965819mHa. The latter starts saving after its first gradient update;
the preserved input seed is better than some early saved outputs. Its actual
best improves that seed slightly and remains invalid.

R28/B64/C32 initializes at1.519059783Ha; Adam12000 reaches
**17.020772619mHa**, below20mHa only before angle quantization. This is a new
floating-point factorization, not native acceptance. Its weight norm is44.68.
Source `grouped_thc_seed.py`/`grouped_thc_fit.py`; output
`candidates/grouped28_b64_c32_float.npz`, logs use the same stem.

`grouped_quantization_probe.py` measures actual decoded signed angle words.
For fixed snapshots, even allowing C64 PSD covariance refits for100cycles gives:

| Snapshot | brot7 error | brot8 error | brot9 error |
|---|---:|---:|---:|
|R24 during L-BFGS|279.057878mHa|171.904133mHa|116.963129mHa|
|R28 during Adam|229.738643mHa|131.728289mHa|not refitted|

Direct rounding errors before those refits are larger. The corresponding
parameter files retain the exact snapshots and quantized vectors. These tests
show that floating-point accuracy alone is insufficient for these dictionaries;
they do not establish global optimality or impossibility for these shapes.
No grouped ANGLE circuits or native submissions were built.
`results/grouped_expanded_shape_heuristic_costs.json` records additional shape
screens, excluding all ANGLES. In particular R28/B64/C32 has heuristicQ358/T7780
atbrot7 butQ629/T6560atbrot8, so higher precision has a material width cost.

No research optimizer remains running. The selected chemistry and archive have
not changed. Other new algorithm experiments are inPHASE_AND_GROUPED_RESEARCH.md.
