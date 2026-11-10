#!/usr/bin/env bash
# Source from anywhere.
ROOT=/media/shared/ExpansionDrive2/QuantumAI
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$ROOT/from_host/kickmix:$ROOT/from_host/verifiers:$ROOT/initial_experiments/4/claude${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1
export PY=/tmp/quantumai-venv.wOrZGC/bin/python
