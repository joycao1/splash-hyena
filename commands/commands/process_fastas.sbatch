#!/bin/bash
#SBATCH --job-name=bkc_faidx
#SBATCH --account=horence
#SBATCH --partition=horence
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/users/%u/bkc_faidx_%A_%a.out
#SBATCH --error=/scratch/users/%u/bkc_faidx_%A_%a.err

set -euo pipefail

MANIFEST="${1:?manifest.txt}"
OUT_ROOT="${2:?OUT_ROOT}"

LOCK="${OUT_ROOT}/.inputs.lock"
LIST="${OUT_ROOT}/dataloader_inputs.txt"

module -q purge
# samtools for faidx; adjust if your cluster uses different module names
module -q load samtools 2>/dev/null || true

# Get my FASTA
idx="${SLURM_ARRAY_TASK_ID:?}"
FA=$(sed -n "$((idx+1))p" "$MANIFEST")
[[ -n "$FA" ]] || { echo "No line $idx in $MANIFEST"; exit 1; }

echo "[$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID] indexing: $FA"

# Index (safe if index exists)
if command -v samtools >/dev/null 2>&1; then
  samtools faidx "$FA"
else
  echo "samtools not found; skipping faidx (ensure *.fai exists)" >&2
fi

# Register path once (atomic)
mkdir -p "$OUT_ROOT"
(
  flock -x 9
  # store absolute path
  realpath "$FA" >> "$LIST"
) 9>"$LOCK"

echo "registered -> $LIST"

