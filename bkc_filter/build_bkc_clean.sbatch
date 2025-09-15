#!/bin/bash
#SBATCH --job-name=bkc_build_clean
#SBATCH --partition=dev
#SBATCH --account=horence
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/users/%u/bkc_build_clean_%j.out
#SBATCH --error=/scratch/users/%u/bkc_build_clean_%j.err

set -euo pipefail
module purge
# Load a newer GCC; tweak if these aren't available on Sherlock
module load gcc/11.3.0 cmake 2>/dev/null || module load gcc/12 2>/dev/null || module load gcc 2>/dev/null || true
export CC="$(which gcc || true)"; export CXX="$(which g++ || true)"
unset LD_LIBRARY_PATH || true

SRC_OAK="/oak/stanford/groups/horence/joycao1/bkc_filter"
BUILD_DIR="/scratch/users/$USER/bkc_build"
BIN_DIR="/scratch/users/$USER/bin"

mkdir -p "$BUILD_DIR" "$BIN_DIR"
rsync -a --delete "$SRC_OAK"/ "$BUILD_DIR"/
cd "$BUILD_DIR"

make clean || true
find . -name '*.o' -delete || true
rm -f libs/mimalloc/mimalloc.o
rm -f libs/zstd/lib/libzstd.a

make -j2 \
  BKC_OUT_BIN_DIR="$BIN_DIR" \
  STATIC_LFLAGS="" \
  LDFLAGS="" \
  CLINK="-lm -std=c++20 -pthread -lz -lzstd" \
  CFLAGS="-fPIC -Wall -O3 -DNDEBUG -std=gnu++20 -pthread \
          -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/mimalloc/include \
          -include stdint.h" \
  CXXFLAGS="-std=gnu++20 -O2 \
            -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/mimalloc/include \
            -include stdint.h" \
  bkc_filter bkc_dump

echo "== Built binaries =="
ls -l "$BIN_DIR"/bkc_filter "$BIN_DIR"/bkc_dump
