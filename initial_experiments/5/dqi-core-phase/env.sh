#!/usr/bin/env bash
CORE_PHASE_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CORE_PHASE_PROJECT=$(cd "$CORE_PHASE_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$CORE_PHASE_ROOT:$CORE_PHASE_PROJECT/initial_experiments/4/codex:$CORE_PHASE_PROJECT/from_host/kickmix:$CORE_PHASE_PROJECT/from_host/verifiers:$CORE_PHASE_PROJECT/initial_experiments/3:$CORE_PHASE_PROJECT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
export TMPDIR="$CORE_PHASE_ROOT/parameters"
