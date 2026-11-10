#!/usr/bin/env bash
QFIT_ROOT=/media/shared/ExpansionDrive2/QuantumAI
export PYTHONPATH="$QFIT_ROOT/initial_experiments/5/chem-qfit:$QFIT_ROOT/initial_experiments/4/codex:$QFIT_ROOT/from_host/kickmix:$QFIT_ROOT/from_host/verifiers:$QFIT_ROOT/initial_experiments/3"
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=3 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export QFIT_PY=/mnt/data01/miniconda3/envs/torch2.13/bin/python
