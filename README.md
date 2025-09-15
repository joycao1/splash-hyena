Optimized Pipeline

This repository contains a high-performance pipeline for processing large single-cell sequencing datasets (e.g. Tabula Sapiens, ERR samples). It includes tools for chunking FASTQs, running the bkc_filter / bkc_dump processing stages, aggregating results, and generating per-cell FASTA outputs ("carrots").

Repository Structure

chunk_fastqs.sh / chunk_fastqs.sbatch
Scripts to split raw FASTQ files into manageable chunks for downstream processing.

chunk_tb_fastqs.sh / chunk_tb_fastqs.sbatch
Specialized chunker for Tabula Sapiens–style lane naming.

pairs.txt / pairs_tabula.txt
Lists of paired FASTQ chunks (R1,R2). Used by the chunking + filter jobs.

run_chunks.sbatch / submit_chunks.sh
Slurm batch scripts to run bkc_filter + bkc_dump in array mode across many chunks.

carrots_ultra.py / carrots_ultra.sbatch
Aggregates per-chunk txt dumps, normalizes, sorts, and converts into FASTA files per sample.

carrots_merge.py / carrots_merge_stream.py
Utilities to merge carrot FASTAs across samples.

count_carrots.sbatch
Helper Slurm job that counts FASTA records (carrots) in each output and reports per-file counts + totals.

commands/
Builds of bkc_filter and bkc_dump executables.

out/
Default output directory for FASTAs, logs, and carrot counts.

Typical Workflow

Chunk FASTQs

sbatch chunk_fastqs.sbatch /path/to/fastqs
# or for Tabula Sapiens
sbatch chunk_tb_fastqs.sbatch /path/to/fastqs


Prepare pairs list
Generated automatically by the chunker (pairs.txt or pairs_tabula.txt).

Run filter/dump on chunks

./submit_chunks.sh pairs.txt
# submits an array of run_chunks.sbatch jobs


Generate FASTAs

NL=$(wc -l < samples.txt)
sbatch --array=1-$NL carrots_ultra.sbatch


Count carrots

sbatch out/count_carrots.sbatch out/fasta carrot_counts

Requirements

Slurm (HPC scheduler)

GNU coreutils (sort, split, awk, etc.)

Python 3.12 (loaded via module load)

bkc_filter / bkc_dump executables (expected in commands/ or PATH)

Optional: pigz, GNU parallel for faster decompression and chunking.

Notes

Adjust --cpus-per-task, --mem, and --time in the sbatch scripts depending on dataset size.

Large FASTA outputs (out/fasta*) are git-ignored; only logs and configs should be committed.

Environment modules (e.g. ml gcc/14.2.0, ml samtools/1.16.1) are loaded inside the batch scripts when needed.
