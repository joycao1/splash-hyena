#!/usr/bin/env bash
set -euo pipefail

# ---------------- user knobs (edit if needed) ----------------
EXEC_FILTER="./commands/bkc_filter"
EXEC_DUMP="./commands/bkc_dump"
ANCHORS="8mer_list_for_single_cell_testing.txt"
# One CSV line "R1,R2" to test scaling on (use a representative pair)
PAIR_LINE="${PAIR_LINE:-}"
# Thread list to test (space separated)
THREADS_LIST="${THREADS_LIST:-4 8 16 32 64 96 128}"
# Output CSV
OUT_CSV="${OUT_CSV:-scaling_results.csv}"

# ---------------- sanity checks ----------------
[[ -x "$EXEC_FILTER" ]] || { echo "ERROR: not executable: $EXEC_FILTER"; exit 1; }
[[ -x "$EXEC_DUMP"   ]] || { echo "ERROR: not executable: $EXEC_DUMP"; exit 1; }
[[ -s "$ANCHORS"     ]] || { echo "ERROR: anchors missing: $ANCHORS"; exit 1; }

if [[ -z "${PAIR_LINE}" ]]; then
  # Try first line of pairs.txt as a fallback
  if [[ -s pairs.txt ]]; then
    PAIR_LINE=$(sed -n '1p' pairs.txt)
  else
    echo "ERROR: set PAIR_LINE='R1,R2' or provide pairs.txt"
    exit 1
  fi
fi

R1=${PAIR_LINE%,*}
R2=${PAIR_LINE#*,}
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: missing input(s): R1=$R1 R2=$R2"; exit 1; }

# ---------------- scratch staging ----------------
SCRATCH_JOB="/scratch/groups/horence/${USER}/sweep_${SLURM_JOB_ID:-$$}"
mkdir -p "${SCRATCH_JOB}"
trap 'rm -rf "${SCRATCH_JOB}"' EXIT

echo "[STAGE] -> ${SCRATCH_JOB}"
rsync -ah "$ANCHORS" "${SCRATCH_JOB}/"
rsync -ah "$R1" "$R2" "${SCRATCH_JOB}/"

ANCHORS_S="${SCRATCH_JOB}/$(basename "$ANCHORS")"
R1_S="${SCRATCH_JOB}/$(basename "$R1")"
R2_S="${SCRATCH_JOB}/$(basename "$R2")"
LIST_S="${SCRATCH_JOB}/pair.csv"
printf '%s,%s\n' "$R1_S" "$R2_S" > "$LIST_S"

# ---------------- utils ----------------
to_seconds() {  # convert h:mm:ss or m:ss to seconds
  awk -F: 'NF==3{print $1*3600+$2*60+$3} NF==2{print $1*60+$2} NF==1{print $1}'
}

measure_run () {
  local nthreads="$1"
  local tag="$2"
  local bkc="${SCRATCH_JOB}/out_${tag}.bkc"
  local txt="${SCRATCH_JOB}/out_${tag}.txt"
  local tfile="${SCRATCH_JOB}/time_${tag}.txt"

  /usr/bin/time -v \
    "$EXEC_FILTER" \
      --mode pair \
      --input_name "$LIST_S" \
      -d "$ANCHORS_S" \
      --cbc_len 16 --umi_len 12 --leader_len 8 --follower_len 31 --gap_len 0 \
      --verbose 1 \
      --n_threads "$nthreads" \
      --output_name "$bkc" \
    1>/dev/null 2>"$tfile"

  /usr/bin/time -v "$EXEC_DUMP" --input_name "$bkc" --output_name "$txt" \
    1>/dev/null 2>>"$tfile"

  # parse /usr/bin/time -v
  local wall=$(grep -m1 "Elapsed (wall clock) time" "$tfile" | awk '{print $8}')
  local rss=$(grep -m1 "Maximum resident set size" "$tfile" | awk '{print $6}')  # KB
  local wall_s=$(echo "$wall" | to_seconds)
  echo "$wall_s" "$rss"
}

# ---------------- warm-up (optional cache warm) ----------------
echo "[WARMUP] 1 run @ 8 threads to warm caches"
measure_run 8 "warmup" >/dev/null

# ---------------- header ----------------
echo "threads,wall_seconds,max_rss_kb" > "$OUT_CSV"

# ---------------- sweep ----------------
for t in $THREADS_LIST; do
  echo "[RUN] threads=${t}"
  # run twice and take the best wall time (reduces noise)
  read -r wall1 rss1 < <(measure_run "$t" "t${t}_a")
  read -r wall2 rss2 < <(measure_run "$t" "t${t}_b")

  # choose best by wall time; keep corresponding rss
  if (( $(printf "%.0f\n" "$wall1") <= $(printf "%.0f\n" "$wall2") )); then
    wall="$wall1"; rss="$rss1"
  else
    wall="$wall2"; rss="$rss2"
  fi
  echo "${t},${wall},${rss}" >> "$OUT_CSV"
done

echo "[DONE] Results -> $OUT_CSV"
echo "Tip: sort -t, -k2,2n $OUT_CSV"

