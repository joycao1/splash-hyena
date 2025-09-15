#!/bin/bash
#SBATCH --job-name=bkc_chunks
#SBATCH --partition=horence
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --open-mode=append

bash script.txt
