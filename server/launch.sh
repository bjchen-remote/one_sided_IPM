#!/usr/bin/env bash
set -euo pipefail

server_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
root="$(dirname -- "$server_dir")"
runner="$root/tools/run_explicit_experiment.py"

usage() {
    cat <<'HELP'
Usage
  bash server/launch.sh smoke [--matlab /path/to/matlab]
  bash server/launch.sh start --case examples/bspline4.json --target 400 --hours 8
  bash server/launch.sh resume /absolute/path/to/run
  bash server/launch.sh status /absolute/path/to/run
Set IPM_RUNS_DIR for the persistent run directory parent and MATLAB_BIN for MATLAB.
Start and smoke accept --threads N. IPM_CHESS_THREADS sets the default (16).
Resume uses the existing run's configured thread count.
HELP
}

case "${1:-}" in
    start|smoke)
        action="$1"
        shift
        stamp="$(date -u +%Y%m%dT%H%M%SZ)"
        run_dir="${IPM_RUNS_DIR:-$server_dir/runs}/ipm_cr_${stamp}_$$_${RANDOM}"
        if [[ "$action" == smoke ]]; then
            options=(start --case "$root/examples/case.json" --threads "${IPM_CHESS_THREADS:-16}" "$@" --target .002 --hours .15)
        else
            options=(start --case "$root/examples/bspline4.json" --threads "${IPM_CHESS_THREADS:-16}" "$@")
        fi
        options+=(--output "$run_dir" --foreground --build-native)
        ;;
    resume|status)
        if [[ $# -ne 2 ]]; then
            usage >&2
            exit 2
        fi
        options=("$1" --output "$2")
        if [[ "$1" == resume ]]; then
            options+=(--foreground)
        fi
        ;;
    -h|--help)
        usage
        exit 0
        ;;
    *)
        usage >&2
        exit 2
        ;;
esac

if [[ -n "${MATLAB_BIN:-}" ]]; then
    options+=(--matlab "$MATLAB_BIN")
elif ! command -v matlab >/dev/null 2>&1 && [[ -x /opt/MATLAB/R2022b/bin/matlab ]]; then
    options+=(--matlab /opt/MATLAB/R2022b/bin/matlab)
fi

exec python3 -B "$runner" "${options[@]}"
