#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source env.sh
"$CENSUS_PY" refine_entry.py --init parameters/transfer_B703.npz --out parameters/repair_B703_bounded.npz --steps 60 --cycles 1 --lambda-cap 320 --seed 5 > results/repair_B703_bounded.log 2>&1
"$CENSUS_PY" refine_entry.py --init parameters/transfer_B703.npz --out parameters/repair_B703_unbounded.npz --steps 60 --cycles 1 --unbounded --seed 5 > results/repair_B703_unbounded.log 2>&1
"$CENSUS_PY" grow.py --start 479 --stop 511 --tag growth511 > results/growth511.log 2>&1
"$CENSUS_PY" grow.py --start 639 --stop 671 --tag growth671 > results/growth671.log 2>&1
"$CENSUS_PY" compiler.py > results/structured.log 2>&1
