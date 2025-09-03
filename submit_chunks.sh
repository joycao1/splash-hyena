#!/usr/bin/env bash
# Usage: THREADS=32 EXEC_FILTER=./bin/bkc_filter EXEC_DUMP=./bin/bkc_dump ./submit_chunks.sh [pairs.txt]
set -euo pipefail
PAIRS="${1:-pairs.txt}"
[[ -s "$PAIRS" ]] || { echo "Error: '$PAIRS' not found or empty"; exit 1; }

NL=$(wc -l < "$PAIRS")
echo "Submitting array with $NL tasks from $PAIRS"
sbatch --array=1-"$NL" run_chunks.sbatch

