DYADIC_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
DYADIC_PROJECT=$(cd "$DYADIC_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$DYADIC_ROOT:$DYADIC_PROJECT/from_host/kickmix:$DYADIC_PROJECT/from_host/verifiers:$DYADIC_PROJECT/initial_experiments/4/codex:$DYADIC_PROJECT/initial_experiments/3:$DYADIC_PROJECT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
