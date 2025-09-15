#!/bin/bash
#SBATCH --job-name=carrots_ultra
#SBATCH --partition=horence
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=05:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

# ---- Config (overridable via env) ----
THREADS="${SLURM_CPUS_PER_TASK:-8}"          # total threads to use for sort
SORT_MEM_GB="${SORT_MEM_GB:-4}"              # per sort buffer size (GNU sort -S)
TXT_DIR="${TXT_DIR:-newtxt}"                 # where the *_R1.part_*.txt live
OUT_DIR="${OUT_DIR:-out/fasta_new}"          # where to write FASTA
SAMPLES_FILE="${SAMPLES_FILE:-samples_tabula.txt}"

module load python/3.12.1
mkdir -p logs "${OUT_DIR}" "${SLURM_TMPDIR:-/tmp}"

# Array line -> sample ID
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p "${SAMPLES_FILE}")
[[ -n "$SAMPLE" ]] || { echo "No sample at line ${SLURM_ARRAY_TASK_ID} in ${SAMPLES_FILE}"; exit 1; }

TMP="${SLURM_TMPDIR:-/tmp}/carrots_${SAMPLE}"
mkdir -p "$TMP"

echo "[INFO] $(date) sample=$SAMPLE"
RAW_N=$(cat "${TXT_DIR}/${SAMPLE}_R1.part_"*.txt | wc -l || true)
echo "[INFO] raw lines: $RAW_N"

# 1) Normalize to 4 cols: cbc \t anchor \t target \t count
# (accept both 5-col (bkc_dump) and 4-col variants)
cat "${TXT_DIR}/${SAMPLE}_R1.part_"*.txt \
| awk '{
    if (NF==5) {printf "%s\t%s\t%s\t%s\n",$2,$3,$4,$5}
    else if (NF==4) {printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}
}' > "$TMP/norm.tsv"

NORM_N=$(wc -l < "$TMP/norm.tsv" || echo 0)
echo "[INFO] normalized lines: $NORM_N"

# 2) Streaming aggregation on disk
#    IMPORTANT: use double quotes so variables expand (fixes your error)
LC_ALL=C sort --parallel="$THREADS" -T "$TMP" -S "${SORT_MEM_GB}G" -k1,1 -k2,2 -k3,3 "$TMP/norm.tsv" \
| awk -F'\t' 'BEGIN{OFS="\t"}
{
  key = $1 OFS $2 OFS $3
  if (key == prev) sum += $4
  else {
    if (NR > 1) print prev, sum
    prev = key; sum = $4
  }
}
END { if (NR > 0) print prev, sum }' \
| LC_ALL=C sort --parallel="$THREADS" -T "$TMP" -S "${SORT_MEM_GB}G" -k1,1 -k2,2 -k4,4nr \
> "$TMP/agg_sorted.tsv"

AGG_N=$(wc -l < "$TMP/agg_sorted.tsv" || echo 0)
echo "[INFO] unique (cbc,anchor,target): $AGG_N"

# 3) Per-CBC stream to FASTA
srun --cpu-bind=cores ./carrots_ultra.py "$TMP/agg_sorted.tsv" "${OUT_DIR}/carrots_${SAMPLE}.fasta"

echo "[DONE] $(date) $SAMPLE -> ${OUT_DIR}/carrots_${SAMPLE}.fasta"

