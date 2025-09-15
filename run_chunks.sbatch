#!/bin/bash
#SBATCH --job-name=bkc_chunks
#SBATCH --partition=horence
#SBATCH --cpus-per-task=32              # ← use the knee from your sweep
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --open-mode=append
# Optional: set an array on submit:  --array=1-N

set -euo pipefail

# --------- Config ----------
ANCHORS="${ANCHORS:-8mer_list_for_single_cell_testing.txt}"
PAIRS_ALL="${PAIRS_ALL:-pairs_tabula.txt}"

# Let caller override thread count; default to the Slurm allocation
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK}}"

# Auto-detect binaries unless provided
EXEC_FILTER="${EXEC_FILTER:-}"
EXEC_DUMP="${EXEC_DUMP:-}"

if [[ -z "${EXEC_FILTER}" ]]; then
  if   [[ -x ./commands/bkc_filter ]]; then EXEC_FILTER=./commands/bkc_filter
  elif [[ -x ./bkc_filter/bin/bkc_filter ]]; then EXEC_FILTER=./bkc_filter/bin/bkc_filter
  elif command -v bkc_filter >/dev/null 2>&1; then EXEC_FILTER=$(command -v bkc_filter)
  else echo "ERROR: bkc_filter not found (try setting EXEC_FILTER=...)"; exit 1
  fi
fi

if [[ -z "${EXEC_DUMP}" ]]; then
  if   [[ -x ./commands/bkc_dump ]]; then EXEC_DUMP=./commands/bkc_dump
  elif [[ -x ./bkc_filter/bin/bkc_dump      ]]; then EXEC_DUMP=./bkc_filter/bin/bkc_dump
  elif command -v bkc_dump >/dev/null 2>&1; then EXEC_DUMP=$(command -v bkc_dump)
  else echo "ERROR: bkc_dump not found (try setting EXEC_DUMP=...)"; exit 1
  fi
fi

# Per-task mini list file
TMP_LIST="$(mktemp -p "${SLURM_TMPDIR:-/tmp}" pairs.XXXX.txt)"
# --------- Env / modules ----------
[[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
module load python
ml gcc/14.2.0

# Avoid srun mem conflicts with allocation
unset SLURM_MEM_PER_CPU SLURM_MEM_PER_GPU
SRUN_MEM_OPT="--mem=0"

# --------- Prep ----------
mkdir -p logs tabuladump tabulatxt "${SLURM_TMPDIR:-/tmp}"

[[ -s "$ANCHORS"   ]] || { echo "ERROR: anchor list not found/empty: $ANCHORS"; exit 1; }
[[ -s "$PAIRS_ALL" ]] || { echo "ERROR: pairs list not found/empty: $PAIRS_ALL"; exit 1; }
[[ -n "${SLURM_ARRAY_TASK_ID-}" ]] || { echo "ERROR: run as an array (--array)"; exit 1; }

PAIR_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PAIRS_ALL" || true)
[[ -n "$PAIR_LINE" ]] || { echo "ERROR: no line ${SLURM_ARRAY_TASK_ID} in $PAIRS_ALL"; exit 1; }

R1=${PAIR_LINE%,*}
R2=${PAIR_LINE#*,}
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: missing input(s): R1=$R1 R2=$R2"; exit 1; }

BASEFILE=$(basename "$R1")
BASEFILE=${BASEFILE%.fastq.gz}; BASEFILE=${BASEFILE%.fq.gz}; BASEFILE=${BASEFILE%.fastq}; BASEFILE=${BASEFILE%.fq}
SAMPLE=$(basename "$(dirname "$R1")")
OUTBASE="${SAMPLE}_${BASEFILE}"

OUT_BKC="newdump/${OUTBASE}.bkc"
OUT_TXT="newtxt/${OUTBASE}.txt"

# Skip if done
if [[ -s "$OUT_BKC" && -s "$OUT_TXT" ]]; then
  echo "[SKIP] Exists: ${OUTBASE}"
  exit 0
fi

# Per-task mini list file
TMP_LIST="$(mktemp -p "${SLURM_TMPDIR:-/tmp}" pairs.XXXX.txt)"
printf '%s,%s\n' "$R1" "$R2" > "$TMP_LIST"

echo "[INFO] JOB=${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} THREADS=${THREADS}"
echo "[INFO] R1=$R1"
echo "[INFO] R2=$R2"
echo "[INFO] OUT_BKC=$OUT_BKC"
echo "[INFO] OUT_TXT=$OUT_TXT"
echo "[INFO] EXEC_FILTER=$EXEC_FILTER"
echo "[INFO] EXEC_DUMP=$EXEC_DUMP"

# --------- Run filter ----------
/usr/bin/time -v srun $SRUN_MEM_OPT --cpu-bind=cores \
  "$EXEC_FILTER" \
    --mode pair \
    --input_name "$TMP_LIST" \
    -d "$ANCHORS" \
    --cbc_len 16 \
    --umi_len 12 \
    --leader_len 8 \
    --follower_len 31 \
    --gap_len 0 \
    --verbose 1 \
    --n_threads "$THREADS" \
    --output_name "$OUT_BKC"

# --------- Dump ----------
/usr/bin/time -v "$EXEC_DUMP" \
  --input_name "$OUT_BKC" \
  --output_name "$OUT_TXT"

echo "[DONE] ${OUTBASE}"

