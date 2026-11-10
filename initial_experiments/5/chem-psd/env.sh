#!/usr/bin/env bash
PSD_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PSD_ROOT=$(cd "$PSD_DIR/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$PSD_DIR:$PSD_ROOT/from_host/kickmix:$PSD_ROOT/from_host/verifiers:$PSD_ROOT/initial_experiments/4/codex:$PSD_ROOT/initial_experiments/5/chem-blocks${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export MPLCONFIGDIR="$PSD_DIR/.matplotlib"
export PSD_PY=/mnt/data01/miniconda3/envs/torch2.13/bin/python
