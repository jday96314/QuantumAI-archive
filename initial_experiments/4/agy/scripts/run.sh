#!/usr/bin/env bash
# Usage: scripts/run.sh script.py [args...]  -- runs with the kickmix/verifier CPU env.
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib
export PYTHONPATH="/media/shared/ExpansionDrive2/QuantumAI/from_host/kickmix:/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers:/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/dqi:/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/scripts:/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy"
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
cd /media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy
exec /tmp/quantumai-venv.wOrZGC/bin/python "$@"
