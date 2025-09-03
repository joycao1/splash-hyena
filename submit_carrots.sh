#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs out/fasta

# Prefer bkctxt/*.txt (from your bkc_dump stage); fall back to out/tsv/*.tsv
if ls -1 bkctxt/*.txt >/dev/null 2>&1; then
  ls -1 bkctxt/*.txt | sort > tsvs.txt
elif ls -1 out/tsv/*.tsv >/dev/null 2>&1; then
  ls -1 out/tsv/*.tsv | sort > tsvs.txt
else
  echo "No TSV-like inputs found (looked in bkctxt/*.txt and out/tsv/*.tsv)"; exit 1
fi

NL=$(wc -l < tsvs.txt)
echo "Submitting carrots array for $NL files"
sbatch --array=1-"$NL" carrots.sbatch

