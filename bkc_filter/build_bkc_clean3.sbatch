#!/bin/bash
#SBATCH --job-name=bkc_build_clean3
#SBATCH --partition=dev
#SBATCH --account=horence
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/users/%u/bkc_build_clean3_%j.out
#SBATCH --error=/scratch/users/%u/bkc_build_clean3_%j.err

set -euo pipefail
module purge
# Use a modern GCC (adjust if Sherlock has a different version)
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

# Hard clean + remove any bad zstd symlink
make clean || true
find . -name '*.o' -delete || true
rm -f libs/mimalloc/mimalloc.o
rm -f libs/zstd/lib/libzstd.a

# 1) Build vendored libs first (no custom C/CXX flags!)
make -j2 libs/zlib-ng/build-g++/zlib-ng/libz.a libs/zstd/lib/libzstd.a

# 2) Now build our two tools with safe flags:
#    - no global CFLAGS injection (avoids breaking zstd asm)
#    - only add -include cstdint to C++.
make -j2 \
  BKC_OUT_BIN_DIR="$BIN_DIR" \
  STATIC_LFLAGS="" \
  LDFLAGS="" \
  CLINK="-lm -std=c++20 -pthread" \
  CXXFLAGS="-std=gnu++20 -O2 -pthread \
            -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/zstd/lib -Ilibs/mimalloc/include \
            -include cstdint" \
  bkc_filter bkc_dump

echo "== Built binaries =="
ls -l "$BIN_DIR"/bkc_filter "$BIN_DIR"/bkc_dump
echo "== ldd check =="
ldd "$BIN_DIR"/bkc_filter || true
ldd "$BIN_DIR"/bkc_dump   || true
