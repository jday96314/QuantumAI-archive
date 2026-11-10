#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source env.sh
# Extend beyond the first draft; preserve the 60-step endpoint separately.
"$CENSUS_PY" refine_entry.py --init parameters/repair_B703_unbounded.npz --out parameters/repair_B703_extended.npz --steps 180 --cycles 1 --unbounded --seed 5 > results/repair_B703_extended.log 2>&1
"$CENSUS_PY" smaller_seed.py > results/smaller_seed.log 2>&1
"$CENSUS_PY" grow.py --start 479 --stop 511 --tag growth511_sequential > results/growth511_sequential.log 2>&1
"$CENSUS_PY" check_components.py > results/structured_components.log 2>&1
