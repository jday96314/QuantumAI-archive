#!/usr/bin/env bash
# Float block fit -> 8-bit integer descent -> offsets -> constrained lambda refit
# -> one-body words -> circuits -> host verification.
# Usage: pipeline.sh FIT.npz TAG COLUMNS [coordinate_steps] [refit_steps]
set -euo pipefail
fit=$1; tag=$2; cols=$3; csteps=${4:-1500}; rsteps=${5:-8000}
here=$(cd "$(dirname "$0")" && pwd); cd "$here"
source ../env.sh
export PYTHONPATH=$here:$PYTHONPATH
PYT=/mnt/data01/miniconda3/envs/torch2.13/bin/python
$PYT block_coordinate.py --init "$fit" --out "q_$tag.npz" --columns "$cols" --steps "$csteps" --moves 64 > "q_$tag.log" 2>&1
$PYT optimize_offsets.py --fit "q_$tag.npz" --out "o_$tag.npz" > "o_$tag.log" 2>&1
$PYT constrained_block.py --fit "o_$tag.npz" --out "c_$tag.npz" --steps "$rsteps" --lr 2e-7 > "cr_$tag.log" 2>&1
$PYT quantize_onebody.py --fit "c_$tag.npz" --out "ob_$tag.npz" --brot 8 --target .012 > "ob_$tag.log" 2>&1
$PY build_global_chem.py --fit "ob_$tag.npz" --out "sub_$tag" > "b_$tag.log" 2>&1
cd /media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers
$PY dfthcsos/verify.py "$here/sub_$tag" > "$here/v_$tag.json" 2> "$here/v_$tag.err" || true
tail -c 600 "$here/v_$tag.json"
