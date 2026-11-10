#!/usr/bin/env bash
set -euo pipefail
cd /media/shared/ExpansionDrive2/QuantumAI/initial_experiments/6/agent14
source ../../5/chem-audit/env.sh
export PYTHONPATH="$PWD:$PYTHONPATH" TMPDIR="$PWD/tmp" XDG_CACHE_HOME="$PWD/cache" CUDA_CACHE_PATH="$PWD/cache/cuda" TORCH_HOME="$PWD/cache/torch" TORCH_EXTENSIONS_DIR="$PWD/cache/torch_extensions" MPLCONFIGDIR="$PWD/cache/matplotlib"
/mnt/data01/miniconda3/envs/torch2.13/bin/python projected_ste_fit.py --init parameters/newton_B671.npz --tag projected_ste_B671 --steps 6000 --lr .01 > results/projected_ste_B671.log 2>&1
