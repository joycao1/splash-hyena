#!/usr/bin/env bash
# wrap.sh — load biology modules and index FASTA files with samtools faidx
# Usage:
#   ./wrap.sh [/path/to/fasta_dir]

set -euo pipefail

# Optional argument for target directory; defaults to current directory
TARGET_DIR="${1:-$PWD}"

# Load modules (Sherlock style)
ml biology
ml samtools/1.16.1

cd "$TARGET_DIR"

# Create .fai indexes for all FASTA files that aren't yet indexed
for f in *.fa *.fasta *.fa.gz *.fasta.gz; do
  [[ -e "$f" ]] || continue         # skip if no match
  if [[ "$f" == *.gz ]]; then
    echo "[SKIP] $f is gzipped—unzip or use samtools faidx directly on uncompressed FASTA"
    continue
  fi
  if [[ -s "$f.fai" ]]; then
    echo "[OK]    Index exists: $f.fai"
  else
    echo "[INFO]  Indexing $f"
    samtools faidx "$f"
  fi
done

echo "All indexing complete."

