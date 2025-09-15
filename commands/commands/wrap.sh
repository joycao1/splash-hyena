#!/bin/bash
#SBATCH --job-name=faidx_wrap
#SBATCH --partition=horence
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
set -euo pipefail

# Always run from submit dir
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Call the wrapper (adjust path if wrap.sh lives elsewhere)
bash ./wrap.sh /scratch/groups/horence/joycao/preprocessing/out/fasta2

