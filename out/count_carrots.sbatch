#!/bin/bash
#SBATCH --job-name=count_carrots
#SBATCH --partition=horence
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=out/logs/%x-%j.out
#SBATCH --error=out/logs/%x-%j.err
#SBATCH --open-mode=append

# Count "carrots" (FASTA records) per file and in total.
# Usage: sbatch out/count_carrots.sbatch [TARGET_DIR] [OUT_LABEL]

set -euo pipefail

# Always run from the submit directory
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Ensure log dir exists for Slurm
mkdir -p out/logs

# Resolve base output dir (use submit dir, not the spool)
OUT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

TARGET_DIR="${1:-${OUT_DIR}/fasta2}"   # or fasta, adjust default if needed
OUT_LABEL="${2:-carrot_counts_fasta2}"
COUNT_DIR="${OUT_DIR}/${OUT_LABEL}"

REPORT_TSV="${COUNT_DIR}/${OUT_LABEL}.tsv"
TOTAL_TXT="${COUNT_DIR}/${OUT_LABEL}_total.txt"

mkdir -p "${COUNT_DIR}"

# Collect candidate files
shopt -s nullglob
FILES=(
  "${TARGET_DIR}"/*.fa
  "${TARGET_DIR}"/*.fasta
  "${TARGET_DIR}"/*.fa.gz
  "${TARGET_DIR}"/*.fasta.gz
)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No FASTA files found under: ${TARGET_DIR}" >&2
  echo "Looked for: *.fa *.fasta *.fa.gz *.fasta.gz" >&2
  exit 1
fi

: > "${REPORT_TSV}"
total=0
for f in "${FILES[@]}"; do
  if [[ "${f}" == *.gz ]]; then
    c=$(zgrep -c '^>' -- "${f}" || true)
  else
    c=$(grep -c '^>' -- "${f}" || true)
  fi
  c=${c:-0}
  printf "%s\t%s\n" "${f}" "${c}" >> "${REPORT_TSV}"
  (( total += c ))
done

echo "${total}" > "${TOTAL_TXT}"

echo "Wrote per-file counts to: ${REPORT_TSV}"
echo "Grand total carrots: ${total}  (also in ${TOTAL_TXT})"
echo
echo "Top 10 by count:"
LC_ALL=C sort -k2,2nr "${REPORT_TSV}" | head -n 10

