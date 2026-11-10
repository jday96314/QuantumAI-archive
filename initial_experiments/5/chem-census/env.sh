#!/usr/bin/env bash
export CENSUS_ROOT=/media/shared/ExpansionDrive2/QuantumAI
export CENSUS_DIR="$CENSUS_ROOT/initial_experiments/5/chem-census"
export CENSUS_PY=/mnt/data01/miniconda3/envs/torch2.13/bin/python
export PYTHONPATH="$CENSUS_DIR:$CENSUS_ROOT/initial_experiments/5/chem-audit:$CENSUS_ROOT/initial_experiments/5/chem-qfit:$CENSUS_ROOT/initial_experiments/4/codex:$CENSUS_ROOT/from_host/kickmix:$CENSUS_ROOT/from_host/verifiers:$CENSUS_ROOT/initial_experiments/3"
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 AUDIT_THREADS=2
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export TMPDIR="$CENSUS_DIR/tmp" CUDA_CACHE_PATH="$CENSUS_DIR/cache/cuda" TORCH_HOME="$CENSUS_DIR/cache/torch"
export XDG_CACHE_HOME="$CENSUS_DIR/cache" MPLCONFIGDIR="$CENSUS_DIR/cache/matplotlib"
