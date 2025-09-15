#!/usr/bin/env bash
# Chunk paired FASTQs (supports Tabula Sapiens lane format).
# Usage:
#   ./chunk_tb_fastqs.sh <INPUT_DIR> [READS_PER_CHUNK=2000000] [OUT_DIR=chunks] [PAIRS_FILE=pairs.txt]
#
# Matches any of:
#   *_R1_*.fastq.gz   (e.g., ..._L001_R1_001.fastq.gz)
#   *_R1.fastq[.gz]
#   *_1.fastq[.gz]
#
# Creates:
#   OUT_DIR/<SAMPLE>/{R1.part_xx.fastq, R2.part_xx.fastq}
#   PAIRS_FILE with lines "<R1chunk>,<R2chunk>"
#
# Notes:
# - SAMPLE is derived so lanes like L001/L002 are preserved (avoids collisions).

set -euo pipefail
shopt -s nullglob

INDIR="${1:-}"
READS_PER_CHUNK="${2:-2000000}"
OUT_DIR="${3:-chunks}"
PAIRS_FILE="${4:-pairs.txt}"

if [[ -z "${INDIR}" || ! -d "${INDIR}" ]]; then
  echo "ERROR: provide INPUT_DIR containing FASTQs" >&2
  echo "Usage: $0 <INPUT_DIR> [READS_PER_CHUNK=2000000] [OUT_DIR=chunks] [PAIRS_FILE=pairs.txt]" >&2
  exit 1
fi

# Validate reads per chunk is a positive integer
if ! [[ "${READS_PER_CHUNK}" =~ ^[0-9]+$ ]] || [[ "${READS_PER_CHUNK}" -le 0 ]]; then
  echo "ERROR: READS_PER_CHUNK must be a positive integer (got: ${READS_PER_CHUNK})" >&2
  exit 1
fi

# lines per FASTQ record chunk (4 lines per read)
LINES=$(( READS_PER_CHUNK * 4 ))

echo "Input dir         : ${INDIR}"
echo "Reads per chunk   : ${READS_PER_CHUNK}  (~${LINES} lines)"
echo "Output chunk root : ${OUT_DIR}"
echo "Pairs list        : ${PAIRS_FILE}"
echo

mkdir -p "${OUT_DIR}"
: > "${PAIRS_FILE}"

# Decompressor selection for .gz
dcmd() {
  if command -v pigz >/dev/null 2>&1; then
    echo "pigz -dc"
  elif command -v zcat >/dev/null 2>&1; then
    echo "zcat"
  else
    echo "gunzip -c"
  fi
}

# Split one FASTQ (gz or not) into OUT_PREFIX_XX.fastq
split_fastq() {
  local in="$1" prefix="$2" lines="$3"
  if [[ "${in}" == *.gz ]]; then
    $(dcmd) "${in}" | split -d -l "${lines}" --additional-suffix=".fastq" - "${prefix}"
  else
    split -d -l "${lines}" --additional-suffix=".fastq" "${in}" "${prefix}"
  fi
}

pairs_found=0

# Loop over candidate R1 files, including Tabula Sapiens lane format
for R1 in \
  "${INDIR}"/*_R1_*.fastq.gz "${INDIR}"/*_R1_*.fq.gz \
  "${INDIR}"/*_R1.fastq     "${INDIR}"/*_R1.fastq.gz \
  "${INDIR}"/*_1.fastq      "${INDIR}"/*_1.fastq.gz
do
  [[ -e "${R1}" ]] || continue

  # Determine mate and SAMPLE name
  if [[ "${R1}" == *"_R1_"* ]]; then
    # e.g. ..._L001_R1_001.fastq.gz → mate: ..._L001_R2_001.fastq.gz
    R2="${R1/_R1_/_R2_}"
    # SAMPLE keeps lane segment (up to _R1_...) to avoid lane collisions
    SAMPLE="$(basename "${R1%%_R1_*}")"
  elif [[ "${R1}" == *_R1.* ]]; then
    R2="${R1/_R1./_R2.}"
    SAMPLE="$(basename "${R1%_R1.*}")"
  else
    # *_1.fastq[.gz]
    R2="${R1/_1./_2.}"
    SAMPLE="$(basename "${R1%_1.*}")"
  fi

  if [[ ! -f "${R2}" ]]; then
    echo "WARN: missing mate for $(basename "${R1}") (expected $(basename "${R2}")); skipping" >&2
    continue
  fi

  pairs_found=$(( pairs_found + 1 ))
  DEST="${OUT_DIR}/${SAMPLE}"
  mkdir -p "${DEST}"

  echo "==> Splitting sample: ${SAMPLE}"
  echo "    R1: $(basename "${R1}")"
  echo "    R2: $(basename "${R2}")"

  split_fastq "${R1}" "${DEST}/R1.part_" "${LINES}"
  split_fastq "${R2}" "${DEST}/R2.part_" "${LINES}"

  # Pair up the new chunks by index
  for C1 in "${DEST}"/R1.part_*.fastq; do
    C2="${C1/R1.part_/R2.part_}"
    if [[ -f "${C2}" ]]; then
      echo "${C1},${C2}" >> "${PAIRS_FILE}"
    else
      echo "WARN: missing chunk mate for $(basename "${C1}")" >&2
    fi
  done
done

if (( pairs_found == 0 )); then
  echo "ERROR: no paired FASTQs found under ${INDIR}" >&2
  echo "       Expected patterns: *_R1_*.fastq.gz, *_R1.fastq[.gz], *_1.fastq[.gz]" >&2
  exit 2
fi

echo
echo "Done. Wrote chunk pairs to: ${PAIRS_FILE}"
wc -l "${PAIRS_FILE}" | awk '{print "Total chunk pairs:", $1}'

