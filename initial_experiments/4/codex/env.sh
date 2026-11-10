#!/usr/bin/env bash
TASK_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$TASK_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$TASK_ROOT:$PROJECT_ROOT/from_host/kickmix:$PROJECT_ROOT/from_host/verifiers:$PROJECT_ROOT/initial_experiments/3:$PROJECT_ROOT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
