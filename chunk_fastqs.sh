#!/usr/bin/env bash
# Chunk all paired FASTQs in a directory into aligned subfiles.
# Usage:
#   ./chunk_fastqs.sh <INPUT_DIR> [READS_PER_CHUNK=2000000] [OUT_DIR=chunks] [PAIRS_FILE=pairs.txt]
#
# Examples:
#   ./chunk_fastqs.sh tdata/dump 500000         # ~0.5M reads per chunk
#   ./chunk_fastqs.sh ./fastq_dumps 250000 chunks_small pairs_small.txt
#
# Notes:
# - Supports: *_1.fastq[.gz], *_2.fastq[.gz], *_R1.fastq[.gz], *_R2.fastq[.gz]
# - Creates:  OUT_DIR/<SAMPLE>/{R1.part_xx.fastq, R2.part_xx.fastq}
# - Outputs:  PAIRS_FILE with lines: "<R1chunk>,<R2chunk>"

set -euo pipefail
shopt -s nullglob

INDIR="${1:-}"
READS_PER_CHUNK="${2:-2000000}"
OUT_DIR="${3:-chunks}"
PAIRS_FILE="${4:-pairs.txt}"

if [[ -z "$INDIR" || ! -d "$INDIR" ]]; then
  echo "ERROR: provide INPUT_DIR containing FASTQs" >&2
  echo "Usage: $0 <INPUT_DIR> [READS_PER_CHUNK=2000000] [OUT_DIR=chunks] [PAIRS_FILE=pairs.txt]" >&2
  exit 1
fi

# lines per FASTQ record chunk (4 lines per read)
LINES=$((READS_PER_CHUNK * 4))

echo "Input dir         : $INDIR"
echo "Reads per chunk   : $READS_PER_CHUNK  (~${LINES} lines)"
echo "Output chunk root : $OUT_DIR"
echo "Pairs list        : $PAIRS_FILE"
echo

mkdir -p "$OUT_DIR"
: > "$PAIRS_FILE"

# helper: split one file (gz or not) into OUT_PREFIX_XX.fastq
split_fastq() {
  local in="$1" prefix="$2" lines="$3"
  if [[ "$in" == *.gz ]]; then
    zcat "$in" | split -d -l "$lines" --additional-suffix=".fastq" - "$prefix"
  else
    split -d -l "$lines" --additional-suffix=".fastq" "$in" "$prefix"
  fi
}

pairs_found=0
for R1 in "$INDIR"/*_1.fastq "$INDIR"/*_1.fastq.gz "$INDIR"/*_R1.fastq "$INDIR"/*_R1.fastq.gz; do
  [[ -e "$R1" ]] || continue

  # determine mate + sample name
  if [[ "$R1" == *_R1.* ]]; then
    R2="${R1/_R1./_R2.}"
    SAMPLE="$(basename "${R1%_R1.*}")"
  else
    R2="${R1/_1./_2.}"
    SAMPLE="$(basename "${R1%_1.*}")"
  fi

  if [[ ! -f "$R2" ]]; then
    echo "WARN: missing mate for $R1 (expected $R2); skipping" >&2
    continue
  fi

  pairs_found=$((pairs_found + 1))
  DEST="$OUT_DIR/$SAMPLE"
  mkdir -p "$DEST"

  echo "==> Splitting sample: $SAMPLE"
  echo "    R1: $(basename "$R1")"
  echo "    R2: $(basename "$R2")"

  # split both reads with identical record counts
  split_fastq "$R1" "$DEST/R1.part_" "$LINES"
  split_fastq "$R2" "$DEST/R2.part_" "$LINES"

  # pair up the new chunks by index
  for C1 in "$DEST"/R1.part_*.fastq; do
    C2="${C1/R1.part_/R2.part_}"
    if [[ -f "$C2" ]]; then
      echo "$C1,$C2" >> "$PAIRS_FILE"
    else
      echo "WARN: missing chunk mate for $(basename "$C1")" >&2
    fi
  done
done

if (( pairs_found == 0 )); then
  echo "ERROR: no *_1/_2 or _R1/_R2 FASTQ pairs found under $INDIR" >&2
  exit 2
fi

echo
echo "Done. Wrote chunk pairs to: $PAIRS_FILE"
wc -l "$PAIRS_FILE" | awk '{print "Total chunk pairs:", $1}'

