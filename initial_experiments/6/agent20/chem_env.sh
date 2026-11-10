#!/usr/bin/env bash
CHEM20_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$CHEM20_ROOT/env.sh"
export EXP_CHEM_PYTHON=/mnt/data01/miniconda3/envs/torch2.13/bin/python
export PYTHONPATH="$CHEM20_ROOT:$CHEM20_ROOT/../../5/chem-audit:$CHEM20_ROOT/../../5/chem-qfit:$PYTHONPATH"
export AUDIT_THREADS=1 MKL_NUM_THREADS=1
export TMPDIR="$CHEM20_ROOT/tmp"
export CUDA_CACHE_PATH="$CHEM20_ROOT/cache/cuda"
export TORCH_HOME="$CHEM20_ROOT/cache/torch"
export XDG_CACHE_HOME="$CHEM20_ROOT/cache"
export MPLCONFIGDIR="$CHEM20_ROOT/cache/matplotlib"
mkdir -p "$TMPDIR" "$CUDA_CACHE_PATH" "$TORCH_HOME" "$MPLCONFIGDIR"
