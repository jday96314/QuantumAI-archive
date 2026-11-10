# Normalization-aware dictionary changes

LATEST checkpoint supersedes the first-cycle numbers below: **9.304538097994794**,
ZIP SHA256`51a94509b0b3b35b191b1ca9f74ecdc89e208a69044743065ccd37e51703d024`.
Selectedchemistry`candidates/norm_cycle2_signed`, fullnativePASS50880,
cost708547388.9553173,lambda48.59417497035283,Q478,T30504. Error1 .011986315443,
error2 .019995304491,Lambda303.5283147432,gap3.9151574724. Canonicaldigest
`ebf8e357a1b693fc8589d5c48a5d6fa10343dfa5944724d8f7a77b8edc2157f7`.
Random-targetchecks2048shotsperall11oracles, exporterrebuild andbotharchiveaudits
PASS. Previous9.2870497 preservedin`revisions/norm_cycle1/`.

Cycle2 was deliberately stopped after79iterations to refit its best-loss
dictionary using the additional physicalerrorroom. The new8000steprefit and
nativepipeline areterminal. The selectedfit changes3505anglewords fromtheJacobi
dictionary, maximumcumulativewordchange118. Frozenfinalfit`parameters/norm_cycle2_final_fit.npz`
SHA256`10688688f0660fd1dc9580ba3541e7db2b6278622cfce48ed18a3c501b11ced6`.
All jobs listed later in this document are nowterminal. Error-only700finished
20.076379mHa; a500-vector latticeprobe finished38.114424mHa, bothinvalid.
See **NORMALIZATION_BOUND.md** for the new global normalization lowerbound and
**DEFERRED_COMPLETION_RESEARCH.md** for the next decoder architecture.
For the selected rebuild commands below substitute `norm_cycle2_final_fit`,
`norm_cycle2_onebody`, `norm_cycle2_plain` and `norm_cycle2_signed` for the
corresponding first-cycle file stems. package_submission.py now selects cycle2.

Current fully verified checkpoint: **9.287049724889583**. The 14.0 goal remains
active. `submission.zip` SHA256 is
`d9409bfab8dc74725bf9963de29d127626c3490732601a815b00af06faae3709`.
The preceding 9.285435820940434 archive and both audits are preserved under
`revisions/signed_gauge/`.

The selected chemistry is `candidates/norm_lattice719_signed`, native result
`results/norm_lattice719_signed_verify.json`: PASS50880, cost712557711.9823189,
lambda_eff48.869214215291805, Q478, T30504. Error1 .011989473875675352,
error2 .0199950128324951, Lambda306.9333489209929, gap3.915394697425654.
Canonical chemistry digest is
`5e47e9d98e4fa67168cdbd460d1430b0a568673067401b3a28c949090f3af352`.
All11 oracles additionally pass2048 random incoming target shots each. The
signed-gauge exporter independently rebuilds all11 KMBs byte-for-byte and all7
NPZ arrays exactly; these are exporter reproducibility checks, not a claim
that the entire numerical optimization was independently replicated.

`package_submission.py` now defaults to this candidate. `audit_submission.py`
independently reopens the final ZIP, binds all14 members to full native results,
recomputes the canonical chemistry digest and overall score, and checks all
default shot counts. Both archive audits pass. Inverse and DQI are unchanged.

## Different objective and representation

The earlier 700-vector lattice search minimizes two-body model error and then
requires a separate normalization refit. `norm_angle_dictionary.py` instead
keeps the existing factor-column orientation. Let P contain packed rank-one
features of the current integer-angle dictionary, and X the frozen physical
factor matrices at the start of a run. For each proposed dictionary use

    w = (P P^T)^-1 P X^T.

The search evaluates the actual reference residual of w^T P, one-body shifts,
normalization and SOS gap. It does not rotate the factor columns into a new
physical eigenbasis. Its objective includes augmented-Lagrangian error/gap
terms; its accepted intermediate states are research fits. Only separately
native-verified candidates enter the archive.

An ordinary unit-word gradient search found no accepted move among192 tested
coordinates. Coordinated LLL/Babai proposals, using the projected feature
curvature, did accept moves. The selected first-cycle fit changes142 angle
words, with maximum cumulative change37 words. It was followed by8000 fixed-
dictionary factor/offset steps at2e-6, then one-body quantization to .012.
The frozen factor input is `parameters/norm_lattice719_final_fit.npz`, SHA256
`cd22bc73a8d3240db5b32fec6af1b9f74070719acd470ecf92ea021375f7b962`.
The earlier gap-projected seed is `parameters/norm_lattice719_gap_native_seed.npz`,
SHA256 `92306f7cc97a070e0fa82ea7f14cd32ea0bc203a7ed1e4c3d293e13e9872561e`.

`--project-gap` uses minimum-norm Newton corrections to the offsets, without
changing error2, to retain a gap margin of3.9155. If X_c are the physical factor
matrices and N_minus the negative-eigenspace projector of h1_prime, the exact
offset derivative is

    d gap / d W_c = W_c + 2 Tr(N_minus X_c).

`results/norm_gap_derivative_check.json` compares it against12 independent
central-difference directions: maximum absolute discrepancy5.07e-8, PASS.

The unprojected first lattice run was deliberately stopped after7 iterations
because it improved loss while going below the native gap floor. Its atomic
best-loss fit was resumed with the gap projection. That run finished after4
iterations, with20 accepted coordinated moves and model-valid lambda48.8902893.
This was a change to enforce the physical gap constraint, not a timeout restart.

## Column generation and exchanges

`augment_norm_dictionary.py` generates atoms from extreme eigenvectors of
matrices dual to each column's L1 representation: solve
`(P P^T) d = sign(w_c)`, and diagonalize the symmetric matrix represented by
`d^T P`. Quantize the extremal vectors and score the actual decoded atoms by
their dual violation and component outside the old span. Sixteen zero-weight
atoms extend the verified719-vector model to735 vectors exactly at initialization.
The generator checks the unchanged old arrays and zero new coefficients.

A16000-step refit reaches model-valid lambda48.6990598. This is not a cost win:
heuristic T rises19900→20266, and a same-style11-oracle family would add176T,
giving estimated cost722693464, above the verified712557712. It remains outside
the archive. This estimate is not a full native result.

`prune_norm_dictionary.py` scores projected row deletions against the full
reference tensor, preserving factor rank and column gauge. Its selected-row
squared-error formula agrees with direct recomputation to below1e-13 in the
recorded probes. Error-only pruning removes all16 new atoms. Forcing all16 new
atoms to remain gives minimum rank350 error20.124914mHa and is invalid even
before normalization. A normalization-aware deletion schedule retains some
new atoms; its12000-step refit is terminal, best model-valid lambda48.9961743,
so it does not beat the selected checkpoint.

`results/separate_spectral_dictionary_probe.json` tests another architecture:
each of350 factor matrices gets its own truncated spectral dictionary (R350,
C1). At B60 the unquantized residual is123.54mHa; B76 recovers the original
two-body model and lowers lambda to47.2332, but needs26676 angle rows and a
larger heuristic width. No circuits were built for this nonwinning probe.
These are results for the tested factor-column orientation, not global bounds.

## Current research jobs

Inspect real processes/logs before deciding whether a job is terminal.

* **PID29022, unified session91305**: second normalization-aware lattice cycle,
  source frozen final fit, all719 rows, maxstep32,512 trials,16 accepted moves
  per iteration,200 requested steps, gap projection, mue .7984,mug .5564.
  Output/log stem `norm_lattice719_cycle2`. At step3 best model-valid lambda
  was48.8135560, below the selected native lambda, but it is not native verified.
* **PID28022, unified session31184**: earlier error-only700-vector continuation,
  output/log stem `lll700_pair_polish_q8_c350`. At step75 best20.0797702mHa,
  still above20mHa. This is not a valid new submission.

Rebuild the selected chemistry from the frozen numerical fit:

```bash
source env.sh
/mnt/data01/miniconda3/envs/torch2.13/bin/python ../../3/quantize_onebody.py \
  --fit parameters/norm_lattice719_final_fit.npz \
  --out parameters/norm_lattice719_final_onebody.npz --brot 8 --steps 3000 --target .012
"$EXP_PYTHON" ../../3/build_global_chem.py \
  --fit parameters/norm_lattice719_final_onebody.npz \
  --out candidates/norm_lattice719_plain --batch 7
"$EXP_PYTHON" build_signed_gauge_chem.py \
  --source candidates/norm_lattice719_plain/submission.npz \
  --out candidates/norm_lattice719_signed --lowbits 8
"$EXP_PYTHON" ../../../from_host/verifiers/dfthcsos/verify.py candidates/norm_lattice719_signed
"$EXP_PYTHON" check_signed_gauge.py --candidate candidates/norm_lattice719_signed
"$EXP_PYTHON" package_submission.py
"$EXP_PYTHON" audit_submission.py
```
