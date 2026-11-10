#!/usr/bin/env bash
set -euo pipefail

if (( $# < 2 )); then
    printf 'Usage: %s {agy|claude|codex} CLI [arguments...]\n' "$0" >&2
    exit 2
fi
agent_name="$1"
shift
agent_executable="$1"
shift

case "$agent_name" in
    agy|claude|codex) ;;
    *) printf 'Unknown agent: %s\n' "$agent_name" >&2; exit 2 ;;
esac

command -v bwrap >/dev/null || { echo 'Install bubblewrap first: sudo apt install bubblewrap' >&2; exit 1; }
resolved_executable="$(command -v "$agent_executable")" || { echo "CLI not found: $agent_executable" >&2; exit 1; }
[[ "$resolved_executable" == /* ]] || { echo "Expected absolute CLI path: $resolved_executable" >&2; exit 1; }

script_directory="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
project_directory="$(cd -- "$script_directory/.." && pwd -P)"
agent_parent="$project_directory/initial_experiments/4"
agent_directory="$agent_parent/$agent_name"

for name in agy claude codex; do
    path="$agent_parent/$name"
    [[ -d "$path" && ! -L "$path" ]] || { echo "Missing directory or symlink: $path" >&2; exit 1; }
done

# Preserve the real HOME, XDG_* settings, credentials, interactive terminal,
# and host runtime sockets. Only hide the other agents' directories.
bubblewrap_arguments=(
    --bind / /
    --unshare-pid
    --proc /proc
    --dev-bind /dev /dev
)
for name in agy claude codex; do
    if [[ "$name" != "$agent_name" ]]; then
        bubblewrap_arguments+=( --perms 000 --tmpfs "$agent_parent/$name" )
    fi
done
bubblewrap_arguments+=(
    --die-with-parent
    --chdir "$agent_directory"
    -- "$resolved_executable" "$@"
)
exec bwrap "${bubblewrap_arguments[@]}"
