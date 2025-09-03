#!/bin/bash
#SBATCH --job-name=carrots_ultra
#SBATCH --partition=horence
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

module load python/3.12.1
mkdir -p logs out/fasta "${SLURM_TMPDIR:-/tmp}"

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p samples.txt)
[[ -n "$SAMPLE" ]] || { echo "No sample at line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

TMP="${SLURM_TMPDIR:-/tmp}/carrots_${SAMPLE}"
mkdir -p "$TMP"

echo "[INFO] $(date) sample=$SAMPLE"
RAW_N=$(cat bkctxt/${SAMPLE}_R1.part_*.txt | wc -l || true)
echo "[INFO] raw lines: $RAW_N"

# 1) Normalize to 4 cols: cbc \t anchor \t target \t count
cat bkctxt/${SAMPLE}_R1.part_*.txt \
| awk '{
    if (NF==5) {printf "%s\t%s\t%s\t%s\n",$2,$3,$4,$5}
    else if (NF==4) {printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$4}
}' > "$TMP/norm.tsv"

NORM_N=$(wc -l < "$TMP/norm.tsv" || echo 0)
echo "[INFO] normalized lines: $NORM_N"

# 2) Streaming aggregation on disk
LC_ALL=C sort --parallel=${SLURM_CPUS_PER_TASK} -T "$TMP" -S 4G -k1,1 -k2,2 -k3,3 "$TMP/norm.tsv" \
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
| LC_ALL=C sort --parallel=${SLURM_CPUS_PER_TASK} -T "$TMP" -S 2G -k1,1 -k2,2 -k4,4nr \
> "$TMP/agg_sorted.tsv"

AGG_N=$(wc -l < "$TMP/agg_sorted.tsv" || echo 0)
echo "[INFO] unique (cbc,anchor,target): $AGG_N"

# 3) Per-CBC stream to FASTA
srun --cpu-bind=cores ./carrots_ultra.py "$TMP/agg_sorted.tsv" "out/fasta/carrots_${SAMPLE}.fasta"

echo "[DONE] $(date) $SAMPLE -> out/fasta/carrots_${SAMPLE}.fasta"

