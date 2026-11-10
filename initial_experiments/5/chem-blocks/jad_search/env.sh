#!/usr/bin/env bash
JAD_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
JAD_PROJECT=$(cd "$JAD_DIR/../../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$JAD_PROJECT/from_host/kickmix:$JAD_PROJECT/from_host/verifiers:$JAD_PROJECT/initial_experiments/4/codex${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1
export JAD_PYTHON=/mnt/data01/miniconda3/envs/torch2.13/bin/python
