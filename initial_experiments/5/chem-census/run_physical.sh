#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
source env.sh
"$CENSUS_PY" prepare_constrained.py > results/continuation703_physical_seed.log 2>&1
"$CENSUS_PY" constrain_entry.py --fit parameters/continuation703_physical_seed.npz --out parameters/continuation703_constrained.npz --steps 6000 --lr 2e-6 --cap .019990 --gap 3.91515 > results/continuation703_constrained.log 2>&1
continuation_fit=parameters/repair_B703_extended.npz
if [ -f parameters/continuation703_constrained.npz ]; then
    continuation_fit=parameters/continuation703_constrained.npz
fi
"$CENSUS_PY" physical_build.py "$continuation_fit" --tag continuation703 > results/continuation703_build.log 2>&1
"$CENSUS_PY" physical_build.py parameters/growth511_sequential_residual_B511.npz --tag growth511 > results/growth511_build.log 2>&1
"$CENSUS_PY" physical_build.py parameters/growth671_residual_B671.npz --tag growth671 > results/growth671_build.log 2>&1
# Salts are declared here before evaluating submissions. No seed override.
flock /tmp/r6-fullverify.lock "$CENSUS_PY" native_record.py candidates/continuation703 --salt round5-chem-census-continuation703 --out continuation703_native.json > results/continuation703_native.log 2>&1
flock /tmp/r6-fullverify.lock "$CENSUS_PY" native_record.py candidates/growth511 --salt round5-chem-census-growth511 --out growth511_native.json > results/growth511_native.log 2>&1
flock /tmp/r6-fullverify.lock "$CENSUS_PY" native_record.py candidates/growth671 --salt round5-chem-census-growth671 --out growth671_native.json > results/growth671_native.log 2>&1
flock /tmp/r6-fullverify.lock "$CENSUS_PY" "$CENSUS_ROOT/from_host/verifiers/dfthcsos/verify.py" candidates/structured --salt round5-chem-census-structured-primary > results/structured_native.json
flock /tmp/r6-fullverify.lock "$CENSUS_PY" "$CENSUS_ROOT/from_host/verifiers/dfthcsos/verify.py" candidates/structured --salt round5-chem-census-structured-independent > results/structured_native_independent.json
