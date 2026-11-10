#!/usr/bin/env bash
set -euo pipefail
script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
exec "$script_directory/run_agent.sh" claude claude --dangerously-skip-permissions "$@"
