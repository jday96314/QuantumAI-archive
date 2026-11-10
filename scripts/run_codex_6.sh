#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 )); then
    printf 'Usage: %s {agent14|agent20} [codex arguments...]\n' "$0" >&2
    exit 2
fi

agent_name="$1"
shift
case "$agent_name" in
    agent14) other_agent_name=agent20 ;;
    agent20) other_agent_name=agent14 ;;
    *) printf 'Unknown agent: %s\n' "$agent_name" >&2; exit 2 ;;
esac

command -v bwrap >/dev/null || { echo 'Install bubblewrap: sudo apt install bubblewrap' >&2; exit 1; }
resolved_executable="$(command -v codex)" || { echo 'codex not found in PATH' >&2; exit 1; }
[[ "$resolved_executable" == /* ]] || { echo "Expected an absolute Codex executable path: $resolved_executable" >&2; exit 1; }

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_directory="$(cd -- "$script_directory/.." && pwd -P)"
agent_parent="$project_directory/initial_experiments/6"
agent_directory="$agent_parent/$agent_name"
other_agent_directory="$agent_parent/$other_agent_name"

for directory in "$agent_directory" "$other_agent_directory"; do
    [[ -d "$directory" && ! -L "$directory" ]] || { echo "Missing directory or symlink: $directory" >&2; exit 1; }
done

# Inherit host filesystem, credentials, networking and terminal. Mask only the
# other agent's working directory with an empty, inaccessible mount.
exec bwrap \
    --bind / / \
    --unshare-pid \
    --proc /proc \
    --dev-bind /dev /dev \
    --perms 000 --tmpfs "$other_agent_directory" \
    --die-with-parent \
    --chdir "$agent_directory" \
    -- "$resolved_executable" --no-daemon --dangerously-bypass-approvals-and-sandbox "$@"
