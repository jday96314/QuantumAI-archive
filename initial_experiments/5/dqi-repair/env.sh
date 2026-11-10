#!/usr/bin/env bash
REPAIR_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$REPAIR_ROOT/../../4/codex/env.sh"
export PYTHONPATH="$REPAIR_ROOT:$PYTHONPATH"
