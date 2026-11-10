#!/usr/bin/env bash
# Verify the delivered artifact bytes without rerunning a stochastic search.
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
source env.sh
"$DISCRETE_PY" "$DISCRETE_ROOT/from_host/verifiers/dfthcsos/verify.py" \
  candidates/final_feasible_B711 --salt round5-chem-discrete-final-feasible-b711 \
  > results/recheck_final_feasible_B711.json
"$DISCRETE_PY" audit_native.py > results/native_digest_audit.log 2>&1
"$DISCRETE_PY" reproduce_native.py > results/native_reproduction.log 2>&1
"$DISCRETE_PY" summarize.py > results/numerical_summary.log 2>&1
