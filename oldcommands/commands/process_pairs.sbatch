#!/bin/bash
#SBATCH --job-name=bkc_pipe
#SBATCH --account=horence
#SBATCH --partition=gpu           # dev/gpu/normal — pick what schedules fast
#SBATCH --gres=gpu:1              # carrots.py uses GPU; set to 0 or remove if CPU
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/%u/bkc_pipe_%A_%a.out
#SBATCH --error=/scratch/users/%u/bkc_pipe_%A_%a.err
# Submit like: sbatch --array=0-(N-1) process_pairs.sbatch manifest.tsv

set -euo pipefail

MANIFEST="${1:?Usage: sbatch --array=0-(N-1) process_pairs.sbatch manifest.tsv}"

# ---- Paths & knobs to edit if needed ----
export PATH="/scratch/users/$USER/bin:$PATH"   # bkc_filter / bkc_dump here

EIGHTMER_DICT="/oak/stanford/groups/horence/joycao1/hyena_preprocess/8mer_list_for_single_cell_testing.txt"

CARROTS="/oak/stanford/groups/horence/joycao1/your_repo/tools/carrots.py"   # <- set this
OUT_ROOT="/oak/stanford/groups/horence/joycao1/hyena_preprocess/outputs"    # <- where to write per-chunk outputs
DATALDR_LIST="${OUT_ROOT}/dataloader_inputs.txt"                            # <- master list for your dataloader
LOCKFILE="${OUT_ROOT}/.dataloader.lock"

# Python / CUDA env (adjust to your cluster)
module purge
module load gcc/12.4.0 2>/dev/null || true
module load cuda 2>/dev/null || true
module load python 2>/dev/null || true
# If you have a venv:
# source /oak/stanford/groups/horence/joycao1/venvs/carrots/bin/activate

# Get my line (array index is 0-based; sed is 1-based)
idx="${SLURM_ARRAY_TASK_ID:?}"
line="$(sed -n "$((idx+1))p" "$MANIFEST")"
[[ -n "$line" ]] || { echo "No line $idx in $MANIFEST"; exit 1; }

IFS=$'\t' read -r R1 R2 RUN_ID CHUNK <<< "$line"

PAIR_TAG="${RUN_ID}_${CHUNK}"
OUT_DIR="${OUT_ROOT}/${RUN_ID}/${CHUNK}"
mkdir -p "$OUT_DIR"

# Make a per-pair fl file (two columns)
FL="${OUT_DIR}/fl.${PAIR_TAG}.txt"
printf "%s\t%s\n" "$R1" "$R2" > "$FL"

# Outputs
BKC="${OUT_DIR}/${PAIR_TAG}.bkc"
TSV="${OUT_DIR}/${PAIR_TAG}.tsv"             # bkc_dump text (TSV with header)
FASTA_DIR="${OUT_DIR}/fasta_out"             # carrots output dir
FASTA="${FASTA_DIR}/carrots.fasta"

echo "== Pair idx=$idx RUN=$RUN_ID CHUNK=$CHUNK =="
echo "R1=$R1"
echo "R2=$R2"
echo "OUT_DIR=$OUT_DIR"

# ---- bkc_filter ----
time bkc_filter \
  --mode pair \
  --input_name "$FL" \
  -d "$EIGHTMER_DICT" \
  --cbc_len 16 \
  --umi_len 12 \
  --leader_len 8 \
  --sample_id "$RUN_ID" \
  --follower_len 31 \
  --gap_len 0 \
  --verbose 2 \
  --output_name "$BKC"

# ---- bkc_dump → TSV (with header) ----
time bkc_dump --input_name "$BKC" --output_name "$TSV"

# ---- carrots.py → fasta_out/carrots.fasta ----
mkdir -p "$FASTA_DIR"
time python "$CARROTS" \
  -i "$TSV" \
  -o "$FASTA_DIR" \
  --output-name carrots.fasta

# ---- register for your dataloader (atomic append) ----
mkdir -p "$OUT_ROOT"
(
  flock -x 9
  echo "$FASTA" >> "$DATALDR_LIST"
) 9>"$LOCKFILE"

echo "Wrote $FASTA and registered in $DATALDR_LIST"

