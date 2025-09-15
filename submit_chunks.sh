#!/usr/bin/env bash
# Usage: THREADS=32 EXEC_FILTER=./bin/bkc_filter EXEC_DUMP=./bin/bkc_dump ./submit_chunks.sh [pairs_tabula.txt]
set -euo pipefail
PAIRS="${1:-pairs_tabula.txt}"
[[ -s "$PAIRS" ]] || { echo "Error: '$PAIRS' not found or empty"; exit 1; }

NL=$(wc -l < "$PAIRS")
echo "Submitting array with $NL tasks from $PAIRS"
sbatch --export=ALL,PAIRS_ALL="$PAIRS" --array=1-"$NL" run_chunks.sbatch

