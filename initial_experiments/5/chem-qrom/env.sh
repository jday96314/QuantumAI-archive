QROM_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
QROM_PROJECT=$(cd "$QROM_ROOT/../../.." && pwd)
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$QROM_ROOT:$QROM_PROJECT/from_host/kickmix:$QROM_PROJECT/from_host/verifiers:$QROM_PROJECT/initial_experiments/4/codex:$QROM_PROJECT/initial_experiments/5/inv-dyadic:$QROM_PROJECT/initial_experiments/3:$QROM_PROJECT/initial_experiments/2${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="$QROM_ROOT/tmp" XDG_CACHE_HOME="$QROM_ROOT/tmp/cache"
export EXP_PYTHON=/tmp/quantumai-venv.wOrZGC/bin/python
