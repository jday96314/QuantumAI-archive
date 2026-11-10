# Multiblock generator and resource search

This isolated subdirectory was created after concurrent changes appeared in
the assigned `5/chem-blocks` directory. It uses only the round-4 Codex-selected
`norm_cycle2_signed` checkpoint, prior-round reference arrays, and host code.
The incumbent remains unchanged and selected. No task-cost win is claimed.

## Exact resource frontier

`resources.py` scans 28,996 constructively costed configurations with brot7/8/9,
chi10, and automatic host `optimize_qubits=True`. It includes ordinary QROMs
and signed-gauge two-entry difference-bank oracles, choosing by the full
`lambda_eff * T * Q` objective. All state-preparation and Givens costs and
the persistent-plus-maximum-temporary width calculation are host-exact.
The full signed symmetric-feature spectrum imposes a necessary PSD rank of
290 for the 20 mHa limit. Each row reports `R*min(B,C)` and the corresponding
optimal-rank residual lower bound; insufficient-rank shapes are excluded from
the reported Pareto frontier. Passing this necessary screen does not establish
physical feasibility.
New ANGLE counts are conservative bounds: every table branch/difference bit
is active and every conditioned Toffoli is charged. Six independent dense
table compilations exactly meet these bounds and stay within 128 qubits.
`results/frontier.csv` contains every scanned row; `results/frontier.json`
contains the Pareto frontier, named shapes, and compilation controls.

Conditional comparison at incumbent lambda_eff = 48.59417497035283:

| Shape | bits/J | host k2 | T bound | Q bound | conditional cost |
|---|---:|---:|---:|---:|---:|
|R1/B719/C350 incumbent|8/7|4|30504|478|708547388.96|
|R2/B351/C175|8/7|4|21863|460|488710645.79|
|R4/B180/C88|7/8|4|17107|443|368266144.19|
|R4/B180/C88|8/7|4|18391|443|395907094.04|
|R4/B224/C88 pivot|8/7|4|21427|443|461263732.48|

These new rows have **untested physical feasibility**, so none is a valid
submission score. R2/B351/C175 is only a conditional 1.45x improvement and
was not fitted. R4/B180 is conditional 1.92x (seven bits) or 1.79x (eight).
R4/B224 at eight bits offers more dictionary freedom while retaining a
conditional 1.54x improvement, motivating the single permitted pivot.

The starting R4/B180 seven-versus-eight-bit k2 example does not hold for the
actual local host. At k2=5 its selection-plus-difference-bank threshold
expression is `ceil(log2(352*ceil(180/32))) + 31*19 = 601`, which fails the
strict `< beta*75` test at both beta7 and beta8. Both choose k2=4. The nearby
R4/B128 case does exhibit a cliff: k2=4 at seven bits and k2=5 at eight;
B129 then returns to k2=4 at eight bits. B256/257 is also explicitly scanned;
for R2/C175 at eight bits, Q rises from 444 to 460 at that boundary.

## Generator derivation and metric checks

Decode the selected integer words first, then form
`P[b,(p,q)] = u[b,p]*u[b,q]*(1 if p=q else sqrt(2))` and
`X[c] = sum_b w[b,c] P[b]`. `unpack(X)` gives the actual symmetric 76x76
generator matrices. The checks compare this with direct `sum_b w*u*u^T`,
the full matrix Frobenius norm, and `P P^T = (u u^T)^2`.

Covariance operations use `P P^T = L L^T`, orthonormal rows
`Q = L^{-1}P`, and physical covariance `H=L^T w w^T L`.
PSD projections are performed in these whitened coordinates. The scripts
check contracted losses against direct subtraction of the full symmetric
feature Gram, retaining the full signed reference spectrum (including its
small negative component). No naive atom-coordinate PSD projection is used.

Orthogonal mixing `Y=O X` also sends offsets to `O offset`. Three deterministic
controls (original, spectral, random orthogonal) preserve both `X^T X` and
`offset^T X` to numerical precision. Approximate commutator distances use
eight fixed Rademacher probes after removing identity components, followed
by balanced column clustering. All groups retain full 76-orbital matrices.
Jacobi joint approximate diagonalization supplies each group's basis.
The final dictionaries may be overcomplete and nonorthogonal: simultaneous
OMP also draws from incumbent atoms and full generator eigenvectors.

The random gauge has the smallest *relative* commutator diagnostic but a
worse physical initializer (1.07309 Ha) than the original gauge (0.75336 Ha).
This shows why the commutator score is an initializer diagnostic rather than
an accuracy or feasibility test. See `results/*generators.json`.

## Bounded fitting and quantization

`relax.py` optimizes unit vectors and whitened weights against the actual
physical two-body Gram residual in float64. Both a mixed-gauge seed and an
original-gauge control are bounded by explicit step and wall-time limits.
Direct matrix subtraction verifies each saved best result.

`quantized_screen.py` uses actual decoded seven/eight-bit words. It tests
nearest signed rounding followed by metric-correct rank-C covariance refits,
and a factorization-aware quantizer with rounded forward passes, differentiable
local weight projection and the full physical tensor objective. Both undergo
further covariance refits. Offset refits preserve the incumbent one-body
correction initially, minimize the attainable SOS gap, and enforce a strict
gap margin where attainable. Actual host `evaluate_words` supplies one-/two-
body errors, Lambda, gap, lambda_eff, and all accuracy margins. These are
native **physical word evaluations**, not ANGLE circuit simulation passes.

Completed results (mHa; each accuracy threshold is 20 mHa):

| Shape/control | continuous error_2 | rounded error_2 | after covariance refit | factor-aware + covariance |
|---|---:|---:|---:|---:|
|R4/B180/C88, mixed gauge, 7-bit screen|65.9524|396.3863|273.7482|273.6894|
|R4/B180/C88, original-gauge continuous control|54.5781|not run|not run|not run|
|R4/B224/C88, single 8-bit pivot|31.8546|262.6675|153.9441|162.9353|

The mixed-gauge fit used 1,801 evaluations (96.5 seconds). The original-gauge
control hit the 250-second limit after 3,307 evaluations; its squared weight
norm grew from about 48.6 to 5955.1 while error plateaued above 54 mHa.
The eight-bit pivot added 44 full-orbital directions per ensemble, preserving
the existing model exactly with zero new weights because the Procrustes
weight projection initially increased error. It used 4,800 evaluations
(177.0 seconds), with a weight ridge of 1e-8 and final squared norm 79.02.
The pivot chooses B224 as the largest nearby scanned four-ensemble dictionary
retaining a conditional >=1.5x resource improvement. No further float search
or second pivot was run.

Actual host physical word evaluations:

| bits/shape/method | error_1 mHa | error_2 mHa | Lambda | E_gap Ha | lambda_eff | T/Q bounds | cost at measured lambda_eff |
|---|---:|---:|---:|---:|---:|---:|---:|
|7/B180/nearest + covariance|129.8269|273.7482|1033.0324|7.3971284|123.4027|17107/443|935195244.86|
|7/B180/factor-aware + covariance|136.4304|273.6894|1035.0055|7.3976363|123.5252|17107/443|936123161.66|
|8/B224/nearest + covariance|61.2822|153.9441|1037.7336|5.7943933|109.5103|21427/443|1039489083.51|
|8/B224/factor-aware + covariance|55.4109|162.9353|984.9990|5.6547576|105.3939|21427/443|1000416322.53|

All four candidates fail **both** accuracy gates. Strict Lambda and E_gap
floors and the valid radicand domain pass, but do not cure failed accuracy.
The finite cost column is an explicitly bounded resource product for an
invalid fit, **not an accepted task score**. The JSONs contain every margin;
`accepted_task_score` is null. No invalid candidate was integrated into a
native chemistry submission, and no new full native chemistry PASS is claimed.
The actual lambda_eff values also invalidate the attractive incumbent-lambda
resource projections for these fits.

`results/oracle_controls.json` records eight separate dense-table XOR controls:
7/8/9 bits, A796 and A972, 64 shots per valid address, two independent simulator
seeds each, arbitrary incoming targets, phase/workspace/selection checks, and
a second application returning the original target. All pass at 120--123
qubits. These are synthetic construction controls, not physical acceptance.
A zero-table padding control also checks that monomial erasure can reserve
workspace when fewer difference bits are active (104Q/37 static Toffolis).

Conclusion: **no replacement; stop after the single pivot**. Preserve the
round-4 native incumbent: lambda_eff 48.59417497035283, T30504, Q478,
cost 708547388.9553173; combined reference 9.3677221531 unchanged.
These results concern the tested initializations, finite budgets, and two
quantization heuristics. They do not establish impossibility of seven-bit or
multiblock factorizations. In particular, the continuous seeds pad two
generator columns with zeros, which remain inactive under the smooth factor
optimization; the subsequent rank-C covariance refits can use the full C88.

## Reproduction

From this subdirectory:

```bash
source ./env.sh
"$JAD_PYTHON" resources.py > results/frontier.log
"$JAD_PYTHON" check_oracle.py > results/oracle_controls.log
"$JAD_PYTHON" generators.py > results/generators.log
"$JAD_PYTHON" generators.py --mix original > results/generators_original.log
"$JAD_PYTHON" relax.py --input candidates/r4_b180_c88_seed.npz --output candidates/r4_b180_c88_float.npz --steps 300 --seconds 600 > results/r4_b180_relax.log
"$JAD_PYTHON" relax.py --input candidates/r4_b180_c88_original_seed.npz --output candidates/r4_b180_c88_original_float.npz --steps 600 --seconds 250 > results/r4_b180_original_relax.log
"$JAD_PYTHON" quantized_screen.py --input candidates/r4_b180_c88_float.npz --bits 7 --cycles 40 --qat-steps 100 > results/r4_b180_quantized.log
"$JAD_PYTHON" grow.py --input candidates/r4_b180_c88_float.npz --B 224 --output candidates/r4_b224_c88_seed.npz > results/r4_b224_grow.log
"$JAD_PYTHON" relax.py --input candidates/r4_b224_c88_seed.npz --output candidates/r4_b224_c88_float.npz --steps 800 --seconds 400 --ridge 1e-8 > results/r4_b224_relax.log
"$JAD_PYTHON" quantized_screen.py --input candidates/r4_b224_c88_float.npz --bits 8 --cycles 60 --qat-steps 150 > results/r4_b224_quantized.log
```

Wall-time limits can change the terminal iteration on another machine;
logs and JSONs record the exact saved evaluations. To reproduce the original
gauge's saved iteration without its machine-dependent cutoff, use
`--steps 551 --seconds 10000`.

Dependencies: NumPy, SciPy, PyTorch, and the bundled kickmix/host verifier.
Read-only imports reuse the round-4 signed-gauge encoder, QROM emitter and
monomial erasure code, with their existing provenance. No external public
algorithm implementation has been copied into this experiment.
