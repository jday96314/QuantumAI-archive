#!/usr/bin/env bash
INV_CENSUS_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
INV_CENSUS_PROJECT=$(cd "$INV_CENSUS_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$INV_CENSUS_ROOT:$INV_CENSUS_PROJECT/from_host/kickmix:$INV_CENSUS_PROJECT/from_host/verifiers:$INV_CENSUS_PROJECT/initial_experiments/3:$INV_CENSUS_PROJECT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export TMPDIR="$INV_CENSUS_ROOT/results"
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
