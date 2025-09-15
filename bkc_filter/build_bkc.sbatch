#!/bin/bash
#SBATCH --job-name=bkc_build
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=bkc_build_%j.out
#SBATCH --error=bkc_build_%j.err

# Load modules
module load gcc cmake 2>/dev/null || true

# Copy to scratch for faster I/O
mkdir -p /scratch/users/$USER/bkc_build
rsync -a --delete /oak/stanford/groups/horence/joycao1/bkc_filter/ /scratch/users/$USER/bkc_build/
cd /scratch/users/$USER/bkc_build

# Clean old objects + remove broken symlink
make clean || true
find . -name '*.o' -delete || true
rm -f libs/mimalloc/mimalloc.o
rm -f libs/zstd/lib/libzstd.a

# Build
make -j4 \
  BKC_OUT_BIN_DIR=/scratch/users/$USER/bin \
  LDFLAGS="" \
  CLINK="-lm -std=c++20 -pthread" \
  CFLAGS="-fPIC -Wall -O3 -DNDEBUG -std=gnu++20 -pthread \
          -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/mimalloc/include \
          -include stdint.h" \
  CXXFLAGS="-std=gnu++20 -O2 \
            -Isrc -Ishared -Ilibs -Ilibs/refresh -Ilibs/mimalloc/include \
            -include stdint.h" \
  bkc_filter bkc_dump

