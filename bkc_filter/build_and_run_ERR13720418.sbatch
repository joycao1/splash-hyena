#!/bin/bash
#SBATCH --job-name=bkc_build_run_ERR13720418
#SBATCH --account=horence
#SBATCH --partition=dev
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:50:00
#SBATCH --output=/scratch/users/%u/bkc_build_run_ERR13720418_%j.out
#SBATCH --error=/scratch/users/%u/bkc_build_run_ERR13720418_%j.err

set -euo pipefail

echo "== [$(date)] Stage 0: toolchain =="
module purge
# Prefer a modern GCC; fall back gracefully
module load gcc/12.4.0 cmake 2>/dev/null || module load gcc/12 cmake 2>/dev/null || module load gcc/11.3.0 cmake 2>/dev/null || module load gcc cmake 2>/dev/null || true
which gcc; gcc --version | head -n1 || true
which g++; g++ --version | head -n1 || true

export CC="$(which gcc || true)"
export CXX="$(which g++ || true)"
unset LD_LIBRARY_PATH || true

# Paths (edit these only if your layout changes)
SRC_OAK="/oak/stanford/groups/horence/joycao1/bkc_filter"
BUILD_DIR="/scratch/users/$USER/bkc_build"
BIN_DIR="/scratch/users/$USER/bin"

FASTQ_DIR="/oak/stanford/groups/horence/joycao1/hyena_preprocess/fastq_data/fastq"
EIGHTMER_DICT="/oak/stanford/groups/horence/joycao1/hyena_preprocess/8mer_list_for_single_cell_testing.txt"
FL_TXT="$FASTQ_DIR/fl.txt"
SAMPLE_ID="ERR13720418"
OUT_BKC="$FASTQ_DIR/${SAMPLE_ID}_pair_1.bkc"
OUT_TXT="$FASTQ_DIR/${SAMPLE_ID}_pair_1.txt"

echo "== [$(date)] Stage 1: rsync source to fast scratch =="
mkdir -p "$BUILD_DIR" "$BIN_DIR"
rsync -a --delete "$SRC_OAK"/ "$BUILD_DIR"/
cd "$BUILD_DIR"

echo "== [$(date)] Stage 2: clean and build vendored libs =="
make clean || true
find . -name '*.o' -delete || true
rm -f libs/mimalloc/mimalloc.o
rm -f libs/zstd/lib/libzstd.a  # clear any stale/dangling symlink

# Build vendored static libs first (no custom flags so their asm/C builds cleanly)
make -j"$(nproc)" libs/zlib-ng/build-g++/zlib-ng/libz.a libs/zstd/lib/libzstd.a

echo "== [$(date)] Stage 3: build bkc_filter & bkc_dump (static libstdc++) =="
# Link statically to libstdc++/libgcc to avoid GLIBCXX mismatches at runtime
make -j"$(nproc)" \
  BKC_OUT_BIN_DIR="$BIN_DIR" \
  STATIC_LFLAGS="" \
  LDFLAGS="" \
  CLINK="-lm -std=c++20 -pthread -static-libstdc++ -static-libgcc" \
  CXXFLAGS="-std=gnu++20 -O2 -pthread \
            -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/zstd/lib -Ilibs/mimalloc/include \
            -include cstdint" \
  bkc_filter bkc_dump

echo "== Built binaries =="
ls -lh "$BIN_DIR"/bkc_filter "$BIN_DIR"/bkc_dump

echo "== ldd check (should NOT show libstdc++.so) =="
ldd "$BIN_DIR"/bkc_filter || true
ldd "$BIN_DIR"/bkc_dump   || true

echo "== [$(date)] Stage 4: run ERR13720418 =="
export PATH="$BIN_DIR:$PATH"
hash -r
which bkc_filter; which bkc_dump

cd "$FASTQ_DIR"

# sanity checks before run
[ -s "$FL_TXT" ] || { echo "Missing input file: $FL_TXT" >&2; exit 2; }
[ -s "$EIGHTMER_DICT" ] || { echo "Missing 8mer dict: $EIGHTMER_DICT" >&2; exit 2; }

echo "Running bkc_filter..."
time bkc_filter \
  --mode pair \
  --input_name "$FL_TXT" \
  -d "$EIGHTMER_DICT" \
  --cbc_len 16 \
  --umi_len 12 \
  --leader_len 8 \
  --sample_id "$SAMPLE_ID" \
  --follower_len 31 \
  --gap_len 0 \
  --verbose 2 \
  --output_name "$OUT_BKC"

echo "Dumping BKC to text..."
time bkc_dump --input_name "$OUT_BKC" --output_name "$OUT_TXT"

echo "== [$(date)] Done. Outputs =="
ls -lh "$OUT_BKC" "$OUT_TXT"
