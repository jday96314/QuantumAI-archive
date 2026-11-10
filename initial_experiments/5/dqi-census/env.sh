#!/usr/bin/env bash
CENSUS_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CENSUS_PROJECT=$(cd "$CENSUS_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$CENSUS_ROOT:$CENSUS_PROJECT/initial_experiments/4/codex:$CENSUS_PROJECT/from_host/kickmix:$CENSUS_PROJECT/from_host/verifiers:$CENSUS_PROJECT/initial_experiments/3:$CENSUS_PROJECT/initial_experiments/2:$CENSUS_PROJECT/initial_experiments/5/dqi-phase${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
