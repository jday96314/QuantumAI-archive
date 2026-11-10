TREE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TREE_ROOT=$(cd "$TREE_DIR/../../.." && pwd)
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1
export LD_LIBRARY_PATH=/mnt/data01/miniconda3/envs/torch2.13/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export PYTHONPATH="$TREE_DIR:$TREE_ROOT/from_host/kickmix:$TREE_ROOT/from_host/verifiers:$TREE_ROOT/initial_experiments/4/codex${PYTHONPATH:+:$PYTHONPATH}"
export TREE_PY=/mnt/data01/miniconda3/envs/torch2.13/bin/python
export TMPDIR="$TREE_DIR/tmp"
export XDG_CACHE_HOME="$TREE_DIR/.cache" CUDA_CACHE_PATH="$TREE_DIR/.cuda-cache"
export MPLCONFIGDIR="$TREE_DIR/.matplotlib"
mkdir -p "$TMPDIR"
