#!/usr/bin/env bash
DISCRETE_ROOT=/media/shared/ExpansionDrive2/QuantumAI
export PYTHONPATH="$DISCRETE_ROOT/initial_experiments/5/chem-discrete:$DISCRETE_ROOT/initial_experiments/4/codex:$DISCRETE_ROOT/from_host/kickmix:$DISCRETE_ROOT/from_host/verifiers:$DISCRETE_ROOT/initial_experiments/3"
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=3 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export DISCRETE_PY=/mnt/data01/miniconda3/envs/torch2.13/bin/python
