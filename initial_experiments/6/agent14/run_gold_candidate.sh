#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source env.sh
"$EXP_PYTHON" build_dp_decoder.py --gold > results/build_gold_decoder.log 2>&1
flock /tmp/r6-fullverify.lock "$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py candidates/gap_gold_n70_m4095.kmb --max-workers 3 > results/gap_gold_verify_default.json 2> results/gap_gold_verify_default.stderr
flock /tmp/r6-fullverify.lock "$EXP_PYTHON" ../../../from_host/verifiers/dqi/verify_submission_dqi.py candidates/gap_gold_n70_m4095.kmb --max-workers 3 --salt agent14-gold-independent-a > results/gap_gold_verify_independent.json 2> results/gap_gold_verify_independent.stderr
"$EXP_PYTHON" ../../5/integration/package_submission.py --output-dir . --dqi candidates/gap_gold_n70_m4095.kmb --dqi-result results/gap_gold_verify_default.json --chemistry ../../5/chem-census/candidates/structured --chemistry-result ../../5/chem-census/results/structured_native.json --previous-summary results/previous_dp_summary.json --reuse-verified-task invsqrt --reuse-verified-task dqi --reuse-verified-task dfthcsos > results/package_gold.log 2>&1
"$EXP_PYTHON" bank_milestones.py > results/bank_gold_milestones.log 2>&1
