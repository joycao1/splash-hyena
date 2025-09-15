# Optimized Pipeline

This repository contains a **high-performance pipeline** for processing **large single-cell sequencing datasets** on the Stanford Sherlock HPC cluster. It includes tools for chunking, filtering, aggregating, and converting sequencing reads into FASTA outputs.

## Repository Structure

- `chunk_fastqs.sh` / `chunk_fastqs.sbatch` – Split raw FASTQ files into manageable chunks.  
- `chunk_tb_fastqs.sh` / `chunk_tb_fastqs.sbatch` – Specialized chunking for Tabula Sapiens lane-format FASTQs.  
- `run_chunks.sbatch` – Slurm batch script for running **bkc_filter** and **bkc_dump** on paired FASTQ chunks.  
- `submit_chunks.sh` – Helper script to submit arrays of chunk-processing jobs.  
- `carrots_ultra.sbatch` / `carrots_ultra.py` – Aggregate filtered results, normalize, sort, and generate FASTA files.  
- `out/count_carrots.sbatch` – Count FASTA records (**“carrots”**) per file and compute grand totals.  
- `samples.txt`, `pairs.txt`, `pairs_tabula.txt` – Input lists for batch jobs.  
- `out/`, `bkcdump/`, `bkctxt/`, `tabuladump/`, `tabulatxt/` – Output and intermediate directories.

## Usage

### 1. Chunk FASTQ Files
```bash
sbatch chunk_fastqs.sbatch /path/to/fastqs
### 2. Prepare Pair and Sample Lists
Ensure pairs.txt or pairs_tabula.txt and samples.txt are correctly generated using sed and sort.

### 3. Run Filtering and Dump
```bash
Copy code
THREADS=32 EXEC_FILTER=./bin/bkc_filter EXEC_DUMP=./bin/bkc_dump ./submit_chunks.sh pairs.txt
### 4. Aggregate to FASTA
```bash
Copy code
NL=$(wc -l < samples.txt)
sbatch --array=1-$NL carrots_ultra.sbatch
### 5. Count Carrots
```bash
Copy code
sbatch out/count_carrots.sbatch out/fasta carrot_counts
### Notes
Partition: Jobs use the horence partition by default.

Modules: Ensure required modules (e.g., python, gcc, samtools) are loaded on Sherlock.

Performance: Increase --cpus-per-task and memory in .sbatch files for very large datasets.

Storage: Outputs are written under out/ (e.g., out/fasta, out/fasta_tabula, carrot_counts_*).

