#!/usr/bin/env bash
CORE_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CORE_PROJECT=$(cd "$CORE_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$CORE_ROOT:$CORE_PROJECT/initial_experiments/4/codex:$CORE_PROJECT/initial_experiments/3:$CORE_PROJECT/initial_experiments/2:$CORE_PROJECT/from_host/kickmix:$CORE_PROJECT/from_host/verifiers${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 PYTHONUNBUFFERED=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
