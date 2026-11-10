#!/usr/bin/env bash
AGENT20_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$AGENT20_ROOT/../../4/codex/env.sh"
export PYTHONPATH="$AGENT20_ROOT:$PYTHONPATH"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
