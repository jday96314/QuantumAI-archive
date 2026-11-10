#!/usr/bin/env bash
AGENT14_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$AGENT14_ROOT/../../4/codex/env.sh"
export PYTHONPATH="$AGENT14_ROOT:$PROJECT_ROOT/initial_experiments/5/dqi-gap:$PROJECT_ROOT/initial_experiments/5/dqi-core/components:$PYTHONPATH"
