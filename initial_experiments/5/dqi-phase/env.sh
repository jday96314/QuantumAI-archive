#!/usr/bin/env bash
PHASE_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PHASE_PROJECT=$(cd "$PHASE_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$PHASE_ROOT:$PHASE_PROJECT/initial_experiments/4/codex:$PHASE_PROJECT/from_host/kickmix:$PHASE_PROJECT/from_host/verifiers:$PHASE_PROJECT/initial_experiments/3:$PHASE_PROJECT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
