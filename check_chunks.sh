#!/usr/bin/env bash
set -euo pipefail

TDATA="/oak/stanford/groups/horence/joycao1/tdata"
CHUNKS="/scratch/groups/horence/joycao/preprocessing/chunks"

printf "%-12s %12s %12s %9s\n" "SAMPLE" "RAW(TDATA)" "CHUNKS" "%COVERED"

for sdir in "$CHUNKS"/ERR*; do
  [[ -d "$sdir" ]] || continue
  s=$(basename "$sdir")

  # Get raw fastq sizes (if they exist)
  if [[ -f "$TDATA/${s}_1.fastq" && -f "$TDATA/${s}_2.fastq" ]]; then
    raw_bytes=$(du -cb "$TDATA/${s}_1.fastq" "$TDATA/${s}_2.fastq" | tail -1 | awk '{print $1}')
  else
    raw_bytes=0
  fi

  # Get chunk sizes
  if ls "$sdir"/*.fastq >/dev/null 2>&1; then
    chunk_bytes=$(du -cb "$sdir"/*.fastq | tail -1 | awk '{print $1}')
  else
    chunk_bytes=0
  fi

  if (( raw_bytes > 0 )); then
    cov=$(( 100 * chunk_bytes / raw_bytes ))
  else
    cov=0
  fi

  # Human readable
  raw_hr=$(numfmt --to=iec --suffix=B "$raw_bytes")
  chk_hr=$(numfmt --to=iec --suffix=B "$chunk_bytes")

  printf "%-12s %12s %12s %8s%%\n" "$s" "$raw_hr" "$chk_hr" "$cov"
done

