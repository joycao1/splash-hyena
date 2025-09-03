#!/bin/bash
#SBATCH --job-name=carrots
#SBATCH --partition=horence
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

# --- env ---
[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
module purge
module load python
ml gcc/14.2.0
# If you use a conda/venv, activate it here:
# source ~/miniconda3/etc/profile.d/conda.sh && conda activate myenv
# or: source /path/to/venv/bin/activate

mkdir -p logs out/fasta

# Pick the Nth TSV line
TSV=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tsvs.txt)
[[ -s "$TSV" ]] || { echo "Missing TSV for task ${SLURM_ARRAY_TASK_ID}: $TSV"; exit 1; }

BASE=$(basename "$TSV"); BASE=${BASE%.tsv}; BASE=${BASE%.txt}  # accept .tsv or .txt

# Run carrots.py from repo root
if [[ -x ./carrots.py ]]; then
  srun --cpu-bind=cores ./carrots.py \
    -i "$TSV" -o out/fasta --output-name "carrots_${BASE}.fasta" \
    > "logs/${BASE}_carrots.log" 2>&1
else
  srun --cpu-bind=cores python ./carrots.py \
    -i "$TSV" -o out/fasta --output-name "carrots_${BASE}.fasta" \
    > "logs/${BASE}_carrots.log" 2>&1
fi

echo "[DONE] carrots ${BASE}"

