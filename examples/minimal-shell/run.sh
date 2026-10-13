#!/usr/bin/env bash
# Copyright (c) 2026 Huang Kai. SPDX-License-Identifier: LGPL-2.1-only
set -euo pipefail

PMTOP_EXAMPLE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PMTOP_ROOT=$(cd -- "$PMTOP_EXAMPLE/../.." && pwd)
if [[ $# -gt 1 || ${1:-} == --help ]]; then
    echo 'Usage: bash examples/minimal-shell/run.sh [new-output-directory]'
    exit $(( $# > 1 ))
fi
for PMTOP_TOOL in cmake icpx icx python3 taskset; do
    if ! command -v "$PMTOP_TOOL" >/dev/null; then
        echo "Missing $PMTOP_TOOL. Load /opt/intel/oneapi/setvars.sh and see the example README." >&2
        exit 1
    fi
done

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OMP_THREAD_LIMIT=1 OMP_DYNAMIC=FALSE MKL_DYNAMIC=FALSE

if [[ $# -eq 1 ]]; then
    if [[ -e $1 ]]; then
        echo 'Output directory already exists; choose a new directory.' >&2
        exit 1
    fi
    mkdir -p -- "$1"
    PMTOP_RUN=$(cd -- "$1" && pwd)
else
    mkdir -p -- "$PMTOP_ROOT/build"
    PMTOP_RUN=$(mktemp -d "$PMTOP_ROOT/build/minimal-shell.XXXXXX")
fi

cmake -S "$PMTOP_ROOT" -B "$PMTOP_ROOT/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx
cmake --build "$PMTOP_ROOT/build" --parallel 2
for PMTOP_INPUT in Config.json Shell_2010_NLIST.lis Shell_2010_ELIST.lis \
                   Shell_2010_ETLIST.lis Shell_2010_DLIST.lis element_library.txt; do
    cp -- "$PMTOP_EXAMPLE/$PMTOP_INPUT" "$PMTOP_RUN/$PMTOP_INPUT"
done
PMTOP_CPU=$(python3 -c 'import os; print(min(os.sched_getaffinity(0)))')
(
    cd -- "$PMTOP_RUN"
    taskset -c "$PMTOP_CPU" "$PMTOP_ROOT/build/src/main/PMTop" 2>&1 | tee stdout.log
)
python3 "$PMTOP_EXAMPLE/validate.py" "$PMTOP_RUN"
printf 'Results: %s\n' "$PMTOP_RUN"
