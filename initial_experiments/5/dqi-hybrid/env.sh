#!/usr/bin/env bash
HYBRID_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
HYBRID_PROJECT=$(cd "$HYBRID_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$HYBRID_ROOT:$HYBRID_PROJECT/initial_experiments/4/codex:$HYBRID_PROJECT/initial_experiments/3:$HYBRID_PROJECT/initial_experiments/2:$HYBRID_PROJECT/from_host/kickmix:$HYBRID_PROJECT/from_host/verifiers${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
