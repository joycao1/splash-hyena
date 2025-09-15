#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-.}"

# Find all R1s ending with ".NNNN.1"; pair to ".NNNN.2"
find "$ROOT" -type f -regextype posix-extended -regex '.*\.[0-9]{4}\.1$' -print \
| sort \
| while IFS= read -r r1; do
  r2="${r1%\.1}.2"
  [[ -f "$r2" ]] || { echo "WARN: missing mate for $r1" >&2; continue; }
  base="$(basename "$r1")"
  run_id="${base%%.*}"
  chunk="$(echo "$base" | awk -F. '{print $(NF-1)}')"
  printf "%s\t%s\t%s\t%s\n" "$r1" "$r2" "$run_id" "$chunk"
done
