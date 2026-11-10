#!/usr/bin/env bash
export AUDIT_ROOT=/media/shared/ExpansionDrive2/QuantumAI
export AUDIT_PY=/mnt/data01/miniconda3/envs/torch2.13/bin/python
export PYTHONPATH="$AUDIT_ROOT/initial_experiments/5/chem-audit:$AUDIT_ROOT/initial_experiments/5/chem-qfit:$AUDIT_ROOT/initial_experiments/4/codex:$AUDIT_ROOT/from_host/kickmix:$AUDIT_ROOT/from_host/verifiers"
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export MPLCONFIGDIR="$AUDIT_ROOT/initial_experiments/5/chem-audit/cache/matplotlib"
export XDG_CACHE_HOME="$AUDIT_ROOT/initial_experiments/5/chem-audit/cache"
export TMPDIR="$AUDIT_ROOT/initial_experiments/5/chem-audit/tmp"
export CUDA_CACHE_PATH="$AUDIT_ROOT/initial_experiments/5/chem-audit/cache/cuda"
export TORCH_HOME="$AUDIT_ROOT/initial_experiments/5/chem-audit/cache/torch"
export TORCH_EXTENSIONS_DIR="$AUDIT_ROOT/initial_experiments/5/chem-audit/cache/torch_extensions"
