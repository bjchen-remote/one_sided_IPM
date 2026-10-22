#!/usr/bin/env bash
set -euo pipefail

# Supply the copied job main.m by absolute path, as in the earlier release.
job_main="${1:?Pass the absolute path of the copied job main.m}"
if [[ "$job_main" != /* || ! -f "$job_main" ]]; then
    echo 'The job entry must be an existing absolute file path.' >&2
    exit 2
fi
case "$job_main" in
    *$'\n'*|*$'\r'*) echo 'The job path must be a single line.' >&2; exit 2 ;;
esac
if [[ -n "${MATLAB_BIN:-}" ]]; then
    matlab_bin="$MATLAB_BIN"
elif command -v matlab >/dev/null 2>&1; then
    matlab_bin="$(command -v matlab)"
elif [[ -x /opt/MATLAB/R2022b/bin/matlab ]]; then
    matlab_bin=/opt/MATLAB/R2022b/bin/matlab
elif [[ -x /opt/MATLAB/R2026a/bin/matlab ]]; then
    matlab_bin=/opt/MATLAB/R2026a/bin/matlab
else
    echo 'Set MATLAB_BIN to the existing server MATLAB executable.' >&2
    exit 127
fi
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export IPM_JOB_MAIN="$job_main"
exec "$matlab_bin" -batch "run(getenv('IPM_JOB_MAIN'))"
