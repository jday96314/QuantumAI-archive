#!/usr/bin/env bash
CHEM_BLOCKS_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
CHEM_PROJECT_ROOT=$(cd "$CHEM_BLOCKS_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$CHEM_BLOCKS_ROOT:$CHEM_PROJECT_ROOT/from_host/kickmix:$CHEM_PROJECT_ROOT/from_host/verifiers:$CHEM_PROJECT_ROOT/initial_experiments/4/codex${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1
export CHEM_PYTHON=/mnt/data01/miniconda3/envs/torch2.13/bin/python
