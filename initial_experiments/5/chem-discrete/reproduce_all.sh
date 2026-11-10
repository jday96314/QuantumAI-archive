#!/usr/bin/env bash
# Long numerical reproduction; writes only this experiment directory.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source env.sh

"$DISCRETE_PY" native_target.py > results/native_target.log 2>&1

"$DISCRETE_PY" reproduce.py > results/reproduction.log 2>&1
"$DISCRETE_PY" "$DISCRETE_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  "$DISCRETE_ROOT/initial_experiments/4/codex/candidates/norm_cycle2_signed" \
  --salt round5-chem-discrete-control > results/incumbent_native.json
"$DISCRETE_PY" rank_frontier.py > results/rank_frontier.log 2>&1
"$DISCRETE_PY" joint_refit.py --init candidates/start_legacy_exact_split_B479.npz \
  --out candidates/joint_B479.npz --steps 300 > results/joint_B479.log 2>&1

"$DISCRETE_PY" admm_refit.py --init candidates/start_legacy_exact_split_B479.npz \
  --out candidates/admm_fullPSD_B479.npz --steps 1800 --rho .01 --full-psd \
  > results/admm_fullPSD_B479.log 2>&1
"$DISCRETE_PY" admm_refit.py --init candidates/admm_fullPSD_B479.npz \
  --out candidates/admm_fullPSD_B479_polish.npz --steps 6000 --rho .003 \
  --full-psd --fixed-rho --warm-dual > results/admm_fullPSD_B479_polish.log 2>&1
"$DISCRETE_PY" psd_certificate.py candidates/admm_fullPSD_B479_polish.npz \
  > results/psd_certificate_polish_B479.log 2>&1
"$DISCRETE_PY" native_pair_control.py > results/native_pair_control.log 2>&1

"$DISCRETE_PY" matched.py > results/matched_driver.log 2>&1
"$DISCRETE_PY" matched_control.py > results/matched_fixed_B479.log 2>&1
"$DISCRETE_PY" matched.py --prefix native_ > results/native_matched_driver.log 2>&1
"$DISCRETE_PY" matched_control.py --prefix native_ > results/native_matched_fixed_B479.log 2>&1
"$DISCRETE_PY" search.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/coupled_B479.npz --mode coupled --steps 200 --seconds 1000 \
  > results/coupled_B479.log 2>&1
"$DISCRETE_PY" search.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/single_B479.npz --mode single --steps 200 --seconds 1000 \
  > results/single_B479.log 2>&1
"$DISCRETE_PY" search.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/thermal_B479.npz --mode thermal --steps 200 --seconds 1000 \
  --temperature .000001 --legacy-temperature > results/thermal_B479.log 2>&1
"$DISCRETE_PY" replace_gpu.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/rankreplace_real_B479.npz --steps 100 --seconds 1000 \
  --trials 4 --temperature .00001 > results/rankreplace_real_B479.log 2>&1
"$DISCRETE_PY" replace.py --init candidates/start_qd_b447_C175.npz \
  --out candidates/basinhop_B447.npz --steps 90 --seconds 1200 --trials 3 \
  --temperature .00001 --quench 4 > results/basinhop_B447_retry.log 2>&1
"$DISCRETE_PY" search.py --init candidates/start_discrete_claude_B479.npz \
  --out candidates/subtree_claude_B479.npz --mode thermal --steps 200 \
  --seconds 1000 --seed 33 --subtree --temperature .000001 --refit-budget 650 \
  --legacy-temperature > results/subtree_claude_B479.log 2>&1

"$DISCRETE_PY" initialize.py --B 479 --seed 0 --rounds 8 > results/init_B479_seed0.log 2>&1
"$DISCRETE_PY" initialize.py --B 447 --seed 1 --rounds 8 > results/init_B447_seed1.log 2>&1
"$DISCRETE_PY" normalize.py --init candidates/refit100_legacy_exact_split_B479.npz \
  --out candidates/norm_B479.npz --steps 4000 --seconds 1500 --lr .00003 \
  > results/norm_B479.log 2>&1
"$DISCRETE_PY" search.py --init candidates/norm_B479_minerror.npz \
  --out candidates/norm_coupled_B479.npz --mode coupled --steps 160 \
  --seconds 1000 --seed 19 --norm-weight .00000001 --refit-budget 650 \
  > results/norm_coupled_B479.log 2>&1
"$DISCRETE_PY" search.py --init candidates/start_pruned_B639.npz \
  --out candidates/coupled_B639.npz --mode coupled --steps 160 --seconds 1000 \
  --seed 7 > results/coupled_B639.log 2>&1
"$DISCRETE_PY" physical.py candidates/coupled_B639.npz --resources \
  > results/b639_control_physical.log 2>&1
"$DISCRETE_PY" search.py --init candidates/evaluated_coupled_B639.npz \
  --out candidates/long_norm_B639.npz --mode coupled --steps 1000 \
  --seconds 1800 --seed 77 --subtree --screen32 --norm-weight .00000001 \
  > results/long_norm_B639.log 2>&1
"$DISCRETE_PY" search.py --init candidates/long_norm_B639_merit.npz \
  --out candidates/native_norm_B639.npz --mode coupled --steps 500 \
  --seconds 450 --seed 78 --subtree --screen32 --norm-weight .00000001 \
  > results/native_norm_B639.log 2>&1
"$DISCRETE_PY" search.py --init candidates/admm_fullPSD_B479_polish.npz \
  --out candidates/psd_geometry_B479_full.npz --mode coupled --steps 250 \
  --seconds 1000 --seed 170 --screen32 \
  --audit-rank-init candidates/start_legacy_exact_split_B479.npz \
  > results/psd_geometry_B479.log 2>&1

"$DISCRETE_PY" normalize.py --init candidates/start_constrained_B711.npz \
  --out candidates/feasible_B711.npz --steps 1800 --seconds 1000 \
  --lr .00003 --jacobi 0 --feasible > results/feasible_B711.log 2>&1
"$DISCRETE_PY" physical.py candidates/feasible_B711.npz --preserve-offset \
  --keep-onebody --resources --build > results/feasible_B711_physical.log 2>&1
"$DISCRETE_PY" "$DISCRETE_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  candidates/feasible_B711 --salt round5-chem-discrete-feasible-b711 \
  > results/feasible_B711_native.json
"$DISCRETE_PY" gauge_escape.py --init candidates/feasible_B711.npz \
  --out candidates/gauge_B711.npz --rounds 10 > results/gauge_B711_retry.log 2>&1
"$DISCRETE_PY" physical.py candidates/gauge_B711.npz --preserve-offset \
  --keep-onebody --resources --build > results/gauge_B711_physical.log 2>&1
"$DISCRETE_PY" "$DISCRETE_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  candidates/gauge_B711 --salt round5-chem-discrete-gauge-b711 \
  > results/gauge_B711_native.json
"$DISCRETE_PY" normalize.py --init candidates/evaluated_gauge_B711.npz \
  --out candidates/final_feasible_B711.npz --steps 200 --seconds 400 \
  --lr .00003 --jacobi 0 --feasible > results/final_feasible_B711.log 2>&1
"$DISCRETE_PY" physical.py candidates/final_feasible_B711.npz --preserve-offset \
  --keep-onebody --resources --build > results/final_feasible_B711_physical.log 2>&1
"$DISCRETE_PY" "$DISCRETE_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  candidates/final_feasible_B711 --salt round5-chem-discrete-final-feasible-b711 \
  > results/final_feasible_B711_native.json

"$DISCRETE_PY" gradient_control.py > results/proposal_gradient_control.log 2>&1
"$DISCRETE_PY" check_oracles.py candidates/gauge_B711/submission.npz \
  candidates/norm_coupled_B479_merit.npz > results/oracle_forward_reverse.log 2>&1
"$DISCRETE_PY" capacity.py candidates/start_legacy_exact_split_B479.npz \
  candidates/start_discrete_claude_B479.npz candidates/start_qd_b447_C175.npz \
  candidates/coupled_B639.npz candidates/norm_coupled_B479.npz \
  candidates/subtree_claude_B479.npz --cpu > results/capacity_final.log 2>&1
"$DISCRETE_PY" summarize.py > results/numerical_summary.log 2>&1
"$DISCRETE_PY" reproduce_native.py > results/native_reproduction.log 2>&1
"$DISCRETE_PY" plot_progress.py > results/plot_progress.log 2>&1
"$DISCRETE_PY" audit_native.py > results/native_digest_audit.log 2>&1
"$DISCRETE_PY" manifest.py > results/manifest.log 2>&1
