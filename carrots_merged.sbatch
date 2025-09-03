#!/bin/bash
#SBATCH --job-name=carrots_merge
#SBATCH --partition=horence
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --open-mode=append

set -euo pipefail

module load python/3.12.1
S=$(sed -n ${SLURM_ARRAY_TASK_ID}p samples.txt)
python3 ./carrots.py --sample "$S" -o out/fasta
