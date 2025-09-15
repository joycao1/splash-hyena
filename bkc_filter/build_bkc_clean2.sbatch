#!/bin/bash
#SBATCH --job-name=bkc_build_clean2
#SBATCH --partition=dev
#SBATCH --account=horence
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/users/%u/bkc_build_clean2_%j.out
#SBATCH --error=/scratch/users/%u/bkc_build_clean2_%j.err

set -euo pipefail
module purge
# Load a modern GCC available on Sherlock; adjust if needed
module load gcc/11.3.0 cmake 2>/dev/null || module load gcc/12 2>/dev/null || module load gcc 2>/dev/null || true
export CC="$(which gcc || true)"
export CXX="$(which g++ || true)"
unset LD_LIBRARY_PATH || true

SRC_OAK="/oak/stanford/groups/horence/joycao1/bkc_filter"
BUILD_DIR="/scratch/users/$USER/bkc_build"
BIN_DIR="/scratch/users/$USER/bin"

mkdir -p "$BUILD_DIR" "$BIN_DIR"
rsync -a --delete "$SRC_OAK"/ "$BUILD_DIR"/
cd "$BUILD_DIR"

# Hard clean
make clean || true
find . -name '*.o' -delete || true
rm -f libs/mimalloc/mimalloc.o
rm -f libs/zstd/lib/libzstd.a  # remove any dangling symlink

# Try to use the system libzstd. Prefer .a; fall back to .so.
LIBZSTD=""
if [ -f /usr/lib64/libzstd.a ]; then
  LIBZSTD=/usr/lib64/libzstd.a
elif [ -f /usr/lib64/libzstd.so ]; then
  LIBZSTD=/usr/lib64/libzstd.so
elif [ -f /usr/lib/x86_64-linux-gnu/libzstd.a ]; then
  LIBZSTD=/usr/lib/x86_64-linux-gnu/libzstd.a
elif [ -f /usr/lib/x86_64-linux-gnu/libzstd.so ]; then
  LIBZSTD=/usr/lib/x86_64-linux-gnu/libzstd.so
else
  echo "Could not find system libzstd; aborting." >&2
  exit 1
fi
echo "Using LIB_ZSTD=$LIBZSTD"

# Build (note: no '-include stdint.h' in CFLAGS; only add '-include cstdint' to CXXFLAGS)
make -j2 \
  BKC_OUT_BIN_DIR="$BIN_DIR" \
  LIB_ZSTD="$LIBZSTD" \
  STATIC_LFLAGS="" \
  LDFLAGS="" \
  CLINK="-lm -std=c++20 -pthread -lz -lzstd" \
  CFLAGS="-fPIC -Wall -O3 -DNDEBUG -std=gnu11 -pthread \
          -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/mimalloc/include" \
  CXXFLAGS="-std=gnu++20 -O2 -pthread \
            -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/mimalloc/include \
            -include cstdint" \
  bkc_filter bkc_dump

echo "== Built binaries =="
ls -l "$BIN_DIR"/bkc_filter "$BIN_DIR"/bkc_dump
echo "== ldd check =="
ldd "$BIN_DIR"/bkc_filter" || true
ldd "$BIN_DIR"/bkc_dump"   || true
