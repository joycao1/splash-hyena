#!/bin/bash
set -euo pipefail
f="$1"; BIN="$2"; DICT="$3"; OUT="$4"; CHUNK_LOG_DIR="$5"
chunk_name="$(basename "$f")"
log="$CHUNK_LOG_DIR/${chunk_name}.log"

# Mirror all stdout/stderr to both Slurm stdout and a per-chunk log
exec > >(tee -a "$log") 2>&1

lines=$(wc -l < "$f" || echo 0)
host=$(hostname)
start_iso=$(date -Is)
echo "[$start_iso] START  $chunk_name  (lines=$lines)  host=$host  pid=$$"

# If empty, skip early but still log it
if [ ! -s "$f" ]; then
  echo "[$(date -Is)] SKIP   $chunk_name (empty chunk)"
  exit 0
fi

# Run bkc_filter; set --verbose 1 for some progress without spamming
/usr/bin/time -v "$BIN" --mode pair \
  --input_name "$f" \
  -d "$DICT" \
  --cbc_len 16 --umi_len 12 --leader_len 8 --follower_len 31 --gap_len 0 \
  --verbose 1 \
  --output_name "$OUT/${chunk_name}.bkc"
rc=$?

end_iso=$(date -Is)
echo "[$end_iso] END    $chunk_name  exit=$rc"
exit $rc
