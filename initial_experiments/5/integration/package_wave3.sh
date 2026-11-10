#!/usr/bin/env bash
set -euo pipefail
WAVE3_TASK_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$WAVE3_TASK_ROOT/env.sh"
cd "$PROJECT_ROOT"

"$EXP_PYTHON" - <<'PY'
import json
from pathlib import Path
path = Path('initial_experiments/5/integration/wave3/results/gap_trimmed_acceptance.json')
assert json.loads(path.read_text())['status'] == 'VERIFIED_IMPROVEMENT'
comparison = path.with_name('new_default_stream_comparison.json')
assert json.loads(comparison.read_text())['status'] == 'PASS'
PY

"$EXP_PYTHON" "$WAVE3_TASK_ROOT/package_submission.py" \
  --output-dir "$WAVE3_TASK_ROOT/wave3" \
  --chemistry "$PROJECT_ROOT/initial_experiments/5/chem-census/candidates/structured" \
  --chemistry-result "$WAVE3_TASK_ROOT/wave3/results/chemistry_verify_default.json" \
  --dqi "$WAVE3_TASK_ROOT/wave3/candidates/gap_trimmed_n70_m4095.kmb" \
  --dqi-result "$WAVE3_TASK_ROOT/wave3/results/gap_trimmed_verify_default.json" \
  --previous-summary "$WAVE3_TASK_ROOT/submission.summary.json" \
  --reuse-verified-task invsqrt --reuse-verified-task dqi --reuse-verified-task dfthcsos

"${WAVE3_METRIC_PYTHON:-/mnt/data01/miniconda3/envs/torch2.13/bin/python}" \
  "$WAVE3_TASK_ROOT/audit_wave3_metric.py"
"$EXP_PYTHON" "$WAVE3_TASK_ROOT/record_wave3_results.py"
