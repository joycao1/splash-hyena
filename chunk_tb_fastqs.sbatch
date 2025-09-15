#!/bin/bash
#SBATCH --job-name=chunk_tb_fastqs
#SBATCH --partition=horence
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

set -euo pipefail

# Always run from the submit directory
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Print env so you know what you got
echo "[$(date -Is)] host=$(hostname) partition=$SLURM_JOB_PARTITION cpus=$SLURM_CPUS_PER_TASK mem=${SLURM_MEM_PER_NODE:-NA}"

# Prefer GNU parallel / pigz if they exist; otherwise the script falls back automatically
if command -v parallel >/dev/null 2>&1; then echo "[info] found: $(command -v parallel)"; else echo "[info] GNU parallel not found; will run serially"; fi
if command -v pigz >/dev/null 2>&1; then echo "[info] found: $(command -v pigz)"; else echo "[info] pigz not found; will use zcat"; fi

# Forward all args to your TB script:
#   ./chunk_tb_fastqs.sh <INPUT_DIR> [READS_PER_CHUNK] [OUT_DIR] [PAIRS_FILE]
bash ./chunk_tb_fastqs.sh "$@"

echo "[$(date -Is)] done"

