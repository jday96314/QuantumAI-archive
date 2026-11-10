#!/usr/bin/env bash
INV_WIDTH_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
INV_PROJECT_ROOT=$(cd "$INV_WIDTH_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$INV_WIDTH_ROOT:$INV_PROJECT_ROOT/from_host/kickmix:$INV_PROJECT_ROOT/from_host/verifiers:$INV_PROJECT_ROOT/initial_experiments/3:$INV_PROJECT_ROOT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
