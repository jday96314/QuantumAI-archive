#!/usr/bin/env bash
# Recheck delivered exact artifacts without repeating numerical searches.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source env.sh
"$TREE_PY" verify_certificates.py > results/independent_dual_audit.log 2>&1
"$TREE_PY" audit.py > results/artifact_audit.log 2>&1
if [[ ${1:-} == --native ]]; then
    "$TREE_PY" "$TREE_ROOT/from_host/verifiers/dfthcsos/verify.py" "$TREE_ROOT/initial_experiments/4/codex/candidates/norm_cycle2_signed" --salt chem-tree-reproduction-control > results/incumbent_native.json 2> results/incumbent_native.stderr
    "$TREE_PY" native_screen.py candidates/grid9_b351_submission candidates/repair_topology_b351_submission candidates/repair_nni_b479_submission > results/reproduced_native.log 2>&1
fi
