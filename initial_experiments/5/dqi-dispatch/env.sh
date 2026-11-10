#!/usr/bin/env bash
DISPATCH_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$DISPATCH_ROOT/../../4/codex/env.sh"
export PYTHONPATH="$DISPATCH_ROOT:$DISPATCH_ROOT/../dqi-repair:$PYTHONPATH"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
