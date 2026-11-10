# Round 5 chemistry multiblock search

**Complete; no replacement.** Two shape families were tested within finite
numerical budgets, including the single eight-bit fallback. Neither meets the
20 mHa physical gates. The 1.5x target was not achieved. The selected chemistry
and combined score remain unchanged. No optimization job from this run remains
running.

Starting checkpoint is the round-4 Codex-selected `norm_cycle2_signed`, read
directly from `../../4/codex/candidates/norm_cycle2_signed/submission.npz`.
Its saved native result passes 50,880 shots: lambda_eff 48.59417497035283,
T 30,504, Q 478, cost 708,547,388.9553173. The combined reference remains
9.3677221531 (inverse 182,956.5272; DQI 1,598,204,229,385).
All incumbent artifacts remain read-only. No replacement has passed yet.

Scope: restructure physical generator ensembles and their shared dictionaries;
do not delete incumbent atoms as the architecture search. Read only host,
shared baseline and prior rounds; write only this directory. No other round-5
agent files are inputs.

## Resource screen

`frontier.py` uses the actual host `SpinFreeData_V1` with chi10 and
`optimize_qubits=True`. It combines persistent widths with the maximum shared
temporary width. Signed-gauge ANGLE costs are constructive worst-case bounds
for `signed_oracle.py`, generalized from the round-4 builder; selected random
table controls are compiled to check those bounds. Physical feasibility is
explicitly untested. Lambda_eff is held at the incumbent only for conditional
resource comparisons; these are not task-cost wins.

Reproduce from this directory:

```bash
source env.sh
"$CHEM_PYTHON" frontier.py > results/frontier.log
```

The host corrects a premise of the exploratory prompt: R4/B180/C88 has k2=4
at both 7 and 8 bits. For b2=19, k2=5 would require 12+31*19=601 temporary
qubits in the host threshold test, exceeding beta*75=600 at eight bits.
R4/B180/C88 has conditional T/Q 17,107/443 (7 bits) and 18,391/443 (8 bits),
giving 1.924x/1.790x only if incumbent lambda and valid physics can be retained.
R2/B351/C175 at eight bits has 21,863/460 and conditional 1.450x.
B256->257 adds a projector-address/word bit and increases workspace; the CSV
and JSON retain both sides. Nine bits frequently crosses to k2=5 and much
larger Q. Ordinary unary QROM is included as an explicitly bounded fallback.

`results/incumbent_native.json` is a fresh full native verification using
independent salt `chem-blocks-baseline-independent-r5`: PASS, 50,880 shots,
identical incumbent task cost and canonical digest. `oracle_controls.py` also
passes 2,048 arbitrary-target shots on each of five synthetic signed-gauge
tables spanning 7/8/9 bits, including a larger address register. Address,
incoming target XOR, zero scratch and absence of phase flips are all checked;
these controls are not physical candidate submissions.

## Physical exploration

Use packed symmetric features with sqrt(2) off-diagonal weights, for which
the Euclidean feature norm equals the full symmetric matrix Frobenius norm.
Gram matrices of rank-one projectors are `(u @ u.T)**2`. Whitening this Gram
is required before any PSD covariance projection. Generator-column rotations
must also rotate identity offsets, preserving both the tensor Gram and the
one-body offset contraction.

`generators.py` derives all 350 full symmetric generator matrices from the
decoded round-4 dictionary. Their packed and full-matrix squared norms both
equal 53.185887057520134. Projector Gram disagreement is 2.11e-15. Direct
reference residual is 19.9953038 mHa, agreeing with native to 7e-7 mHa.
Incumbent, physical-PCA and seeded random orthogonal column gauges are checked;
Gram differences are <=2.34e-15, offset contractions <=4.27e-14, and offset
squared norms <=1.82e-12. All inter-orbital entries remain present.

Balanced commutator clustering and QR-based joint approximate diagonalization
produce per-ensemble bases. JAD is too restrictive: at R4 the chosen random
gauge leaves 32.7%-42.6% relative factor offdiagonal norms. The expanded atom
pools combine JAD vectors, full generator eigenvectors, and the complete
incumbent dictionary; symmetric-metric simultaneous OMP selects nonorthogonal
overcomplete bases. This is a generator/subspace reorganization rather than
an orbital-disjoint split or simple incumbent-row pruning.

Bounded float64 full-tensor optimization uses Cholesky whitening of the
projector Gram, retains negative target modes, and checks contracted loss
against direct tensor subtraction. R4/B180/C88 Adam 2,400 steps reaches
88.134 mHa from the random-column warm start and 58.075 mHa from PCA. Fixed
300-step L-BFGS refinements follow; no seven-bit impossibility is inferred.
The one permitted fallback family is R2/B351/C175, eight-bit targeted.

Final continuous results (all **unquantized**, not submissions):

| Shape | Warm start | Adam 2,400 error_2 (mHa) | L-BFGS 300 error_2 (mHa) |
|---|---|---:|---:|
| R4/B180/C88 | Random orthogonal columns | 88.134 | 63.627 |
| R4/B180/C88 | Physical PCA columns | 58.075 | 54.961 |
| R2/B351/C175 | Physical PCA columns | 29.658 | 27.153 |

The best R4 PCA refinement grows weight squared norm to 1,077.98 while
remaining invalid. The lower-weight Adam snapshots are retained and screened
separately so that this ill-conditioning does not become the sole precision
control. These are local numerical results, not capacity or precision bounds.

## Quantized physical gates

`quantized_screen.py` decodes actual signed-gauge angle words, performs 40
rank-C PSD block-refit cycles in the whitened symmetric feature metric,
tries up to 32 rounds of gradient-ranked integer coordinate proposals against
the **full physical tensor objective**, then performs another 40 block-refit
cycles. Offset refits preserve the associated symmetric contraction in least
squares against the incumbent target. One-body angles are computed for the
resulting Hamiltonian. Native physical evaluation checks both errors,
normalization, gap and domain. An extra integer-move identity control verifies
four accepted updates against independently reconstructed full tensors.

| Shape / snapshot | Bits | error_1 mHa | error_2 mHa | Lambda | E_gap Ha | lambda_eff |
|---|---:|---:|---:|---:|---:|---:|
| R4 PCA Adam | 7 | 133.754 | 365.655 | 1,838.221 | 360.382424 | 1,093.183 |
| R4 PCA L-BFGS | 7 | 135.043 | 353.211 | 1,639.145 | 200.345628 | 785.272 |
| R4 PCA L-BFGS | 8 | 60.116 | 216.362 | 1,559.289 | 96.611548 | 540.330 |
| R2 PCA Adam | 8 | 52.202 | 179.874 | 1,645.821 | 7.452177 | 156.443 |
| R2 PCA L-BFGS | 8 | 55.478 | 172.525 | 1,641.430 | 6.611230 | 147.173 |

Every row fails **both** error gates. All five pass the strict Lambda>49,
E_gap>3.914898 and gap-domain checks. Exact signed margins in Ha, including
every strict-gate margin, are stored in the corresponding
`results/*_quantized.json` and consolidated in `results/final_audit.json`.
For example the best fallback has error margins -35.477842/-152.524994 mHa,
Lambda margin +1,592.429520, and gap margin +2.696332 Ha. Its effective
normalization is also far above the incumbent. Further one-body/normalization
polishing cannot rescue its failed two-body gate, so it was not pursued.

These failures follow decoded quantization, constrained weight refits, integer
angle improvements, and multiple snapshots; they do not establish seven-bit
impossibility. No failed family was exported as an accepted candidate or
assigned a native task cost. The T/Q figures above remain explicitly bounded
resource projections. No new oracle-family/native submission verification was
warranted. The generalized oracle implementation itself has native synthetic
XOR/phase/cleanup validation under the 128-qubit cap.

## Reference precision correction

The cached prior-round training target is not bit-identical to the host:
max entry difference 2.9802324719163664e-8 Ha and Frobenius difference
3.1116404560991105e-7 Ha (0.000311164 mHa). The first quantization preflight
rejected a too-tight 1e-9 assertion. After measuring the discrepancy, all
quantized objectives and final physical metrics use the **fresh native
FCIDUMP and fixed-shift target**, including its negative modes. No host file
was changed. Continuous fits retain the declared cached target; their errors
have at most 0.000311164 mHa target uncertainty, far below their failed margins.
The rank-289/290 errors and the 1-Lipschitz distance bound certify that this
discrepancy does not change the frontier's necessary unrestricted PSD rank290.

## Final selection and artifacts

The only complete valid submission remains the original R1/B719/C350,
8-bit/J7/chi10: lambda_eff **48.59417497035283**, native T **30,504**, native
Q **478**, cost **708,547,388.9553173**, chemistry improvement **1.000x**.
It has one-body margin +8.013684557 mHa, two-body margin +0.004695509 mHa,
Lambda margin +254.528314743, and gap margin +0.259472366 mHa.

`results/final_audit.json` binds its twelve chemistry files to SHA-256 hashes,
checks the fresh native result against the round-4 selected canonical digest,
hashes the research NPZ files/builders/evidence, and records all rejected
screens. It does not replace or repackage the preserved round-4 submission.
`parameters/` contains generator gauges, continuous fits and decoded quantized
research NPZ files. `signed_oracle.py` is the generalized oracle builder.
`results/frontier.csv` contains 86,552 configurations; `frontier.json` includes
19 T/Q Pareto points and detailed B255/256/257 scan evidence (CSV for B255).
All frontier feasibility labels remain untested until physical fitting.

Reused local implementations are explicitly imported from round4:
`angle_gauge_probe.encode`, `build_selectswap_chem.emit_qrom`, and
`monomial_unlookup.erase`. Host resource and physical evaluators are imported
unchanged from `from_host/verifiers`; their Apache-2.0 license is retained at
`../../../from_host/verifiers/LICENSE.txt`. No additional public code or DQI
decoder code is used in this experiment. Runtime: Python3.13, torch2.13.0,
NumPy2.5.2, SciPy1.18.0, and the bundled native kickmix extension.

## Reproduction

Run from `initial_experiments/5/chem-blocks` after `source env.sh`.

```bash
"$CHEM_PYTHON" frontier.py > results/frontier.log
"$CHEM_PYTHON" oracle_controls.py > results/oracle_controls.log
"$CHEM_PYTHON" ../../../from_host/verifiers/dfthcsos/verify.py \
  ../../4/codex/candidates/norm_cycle2_signed \
  --salt chem-blocks-baseline-independent-r5 > results/incumbent_native.json

"$CHEM_PYTHON" generators.py --R 4 --B 180 --C 88 --jad-steps 400 \
  --tag r4b180c88 > results/r4b180c88_generators.log
"$CHEM_PYTHON" generators.py --R 4 --B 180 --C 88 --jad-steps 400 \
  --gauge physical_PCA --tag r4b180c88_pca > results/r4b180c88_pca_generators.log
"$CHEM_PYTHON" generators.py --R 2 --B 351 --C 175 --jad-steps 400 \
  --gauge physical_PCA --tag r2b351c175 > results/r2b351c175_generators.log

for tag in r4b180c88 r4b180c88_pca r2b351c175; do
  "$CHEM_PYTHON" relax.py --input "parameters/${tag}_seed.npz" \
    --output "parameters/${tag}_adam.npz" --steps 2400 --lr .002 --trace \
    > "results/${tag}_adam.log"
  "$CHEM_PYTHON" relax.py --input "parameters/${tag}_adam.npz" \
    --output "parameters/${tag}_lbfgs.npz" --steps 300 --lr 1 \
    --optimizer lbfgs --trace > "results/${tag}_lbfgs.log"
done

"$CHEM_PYTHON" quantized_screen.py --input parameters/r4b180c88_pca_lbfgs.npz \
  --bits 7 8 --cycles 40 --moves 32 > results/r4b180c88_quantized.log
"$CHEM_PYTHON" quantized_screen.py --input parameters/r4b180c88_pca_adam.npz \
  --bits 7 --cycles 40 --moves 32 > results/r4b180c88_adam_quantized.log
for snapshot in adam lbfgs; do
  "$CHEM_PYTHON" quantized_screen.py \
    --input "parameters/r2b351c175_${snapshot}.npz" --bits 8 --cycles 40 \
    --moves 32 > "results/r2b351c175_${snapshot}_quantized.log"
done
"$CHEM_PYTHON" integer_metric_control.py > results/integer_metric_control.log
"$CHEM_PYTHON" final_audit.py > results/final_audit.log
```

Floating-point optimization timings and last-digit trajectories can depend on
GPU scheduling; stored NPZ hashes identify the exact reported snapshots. All
acceptance comparisons here use the stored decoded candidates and native
physical evaluator. No DQI or inverse candidate was modified or retested.


## Isolated JAD search completed

The search under [jad_search/RESULTS.md](jad_search/RESULTS.md) is complete.
It records a 28,996-row host-exact resource scan with explicit ANGLE upper
bounds, metric-correct generator reconstruction, three orthogonal-mixing
controls, full-orbital column clustering/JAD, overcomplete dictionaries, and
bounded relaxations. R4/B180/C88 continuous controls reached 65.95/54.58 mHa;
its best tested seven-bit covariance/quantizer control was 273.69 mHa. The
single eight-bit R4/B224/C88 pivot reached 31.85 mHa continuously and 153.94
mHa in the best tested quantized two-body check. All tested quantized
candidates fail both accuracy gates. No replacement or task-cost win.
Round-4 Codex 708,547,388.9553173 chemistry and combined 9.3677221531 remain
selected. Search stopped after this one pivot; no further optimizer is running.

The local host chooses k2=4 for R4/B180/C88 at **both** seven and eight bits;
the starting seven/eight-bit cliff example does not apply to that shape.
Synthetic signed XOR controls pass 64 shots/address with two independent
simulator seeds, at 120--123 qubits, but are not chemistry submissions.
Native physical-word metrics and every acceptance margin are in
[jad_search/results/FINAL.json](jad_search/results/FINAL.json); none is a full
native submission PASS. Reproduction commands, failed NPZs, builders, logs
and artifact hashes are preserved under `jad_search/`.
