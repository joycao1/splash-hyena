#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Merge multiple bkctxt chunk files for one sample into a single FASTA.

Input files: bkctxt/<SAMPLE>_R1.part_XX.txt
- Each line has either 4 or 5 columns, whitespace-separated (tab/space):
    4-col: cbc  anchor  target  count
    5-col: <ignored>  cbc  anchor  target  count   (first column often 0)
- We GROUP BY (cbc, anchor, target) and SUM counts across all chunks
- Then drop targets that are "universal" across all anchors for a CBC (when CBC has >1 anchors)
- Finally, build one concatenated sequence per CBC:
    for each anchor (in sorted order as read), append anchor once, then all its targets (sorted by count desc)

Output: out/fasta/carrots_<SAMPLE>.fasta
"""

import argparse
from pathlib import Path
from typing import List

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_args():
    ap = argparse.ArgumentParser(
        description="Merge bkctxt chunks for a sample into FASTA (sum counts across chunks)."
    )
    ap.add_argument("-s", "--sample", required=True,
                    help="Sample ID prefix, e.g. ERR13720419 (matches bkctxt/<S>_R1.part_*.txt)")
    ap.add_argument("-i", "--input-dir", default="bkctxt",
                    help="Directory containing chunk .txt files (default: bkctxt)")
    ap.add_argument("-o", "--output-dir", default="out/fasta",
                    help="Directory to write FASTA (default: out/fasta)")
    ap.add_argument("--print-head", type=int, default=0,
                    help="Print first N rows of merged table for sanity (0 disables)")
    return ap.parse_args()


def glob_inputs(input_dir: Path, sample: str) -> List[Path]:
    return sorted((input_dir).glob(f"{sample}_R1.part_*.txt"))


def read_chunk(path: Path) -> pd.DataFrame:
    """
    Read a single chunk file that may be 4 or 5 columns, whitespace-separated.
    Normalize to columns: cbc, anchor, target, count (count=int).
    """
    # Try flexible whitespace sep, no header
    df = pd.read_csv(path, sep=r"\s+", header=None, engine="python")
    if df.shape[1] == 5:
        df.columns = ["_drop", "cbc", "anchor", "target", "count"]
        df = df[["cbc", "anchor", "target", "count"]]
    elif df.shape[1] == 4:
        df.columns = ["cbc", "anchor", "target", "count"]
    else:
        raise ValueError(f"{path}: expected 4 or 5 columns, got {df.shape[1]}")

    # Normalize types
    df["cbc"] = df["cbc"].astype(str)
    df["anchor"] = df["anchor"].astype(str)
    df["target"] = df["target"].astype(str)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)

    # Drop empties just in case
    df = df[(df["anchor"].str.len() > 0) & (df["target"].str.len() > 0)].copy()
    return df


def load_and_merge(files: List[Path]) -> pd.DataFrame:
    if not files:
        raise FileNotFoundError("No input chunk files matched for this sample.")
    frames = [read_chunk(f) for f in files]
    big = pd.concat(frames, ignore_index=True)
    merged = (
        big.groupby(["cbc", "anchor", "target"], as_index=False)["count"]
           .sum()
           .sort_values(["cbc", "anchor", "count"], ascending=[True, True, False])
           .reset_index(drop=True)
    )
    return merged


def filter_universal_targets(df_sorted: pd.DataFrame) -> pd.DataFrame:
    """
    Drop rows where (cbc, target) appears in ALL anchors for that CBC (only when CBC has >1 anchors).
    """
    n_anchors_per_cbc = df_sorted.groupby("cbc")["anchor"].nunique().rename("n_anchors")
    anchors_per_target = (
        df_sorted.groupby(["cbc", "target"])["anchor"]
                 .nunique()
                 .rename("anchors_for_target")
                 .reset_index()
    )
    tmp = df_sorted.merge(anchors_per_target, on=["cbc", "target"], how="left")
    tmp = tmp.merge(n_anchors_per_cbc, on="cbc", how="left")

    keep = ~((tmp["anchors_for_target"] == tmp["n_anchors"]) & (tmp["n_anchors"] > 1))
    return tmp.loc[keep, ["cbc", "anchor", "target", "count"]]


def build_sequences(df_sorted: pd.DataFrame) -> dict:
    """
    For each CBC:
      for each anchor (current order), append anchor once, then all targets for that anchor (current order).
    """
    cbc_to_seq = {}
    for cbc, g in df_sorted.groupby("cbc", sort=False):
        parts = []
        for anchor, a in g.groupby("anchor", sort=False):
            parts.append(anchor)
            parts.extend(a["target"].tolist())
        cbc_to_seq[cbc] = "".join(parts)
    return cbc_to_seq


def write_fasta(cbc_to_seq: dict, out_path: Path) -> int:
    def sanitize(s: str) -> str:
        allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
        return "".join(ch if ch in allowed else "_" for ch in s)
    recs = [SeqRecord(Seq(seq), id=sanitize(cbc), description="") for cbc, seq in cbc_to_seq.items()]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return SeqIO.write(recs, str(out_path), "fasta")


def main():
    args = parse_args()
    files = glob_inputs(Path(args.input_dir), args.sample)
    merged = load_and_merge(files)

    if args.print_head > 0:
        print(merged.head(args.print_head).to_string(index=False))

    filtered = filter_universal_targets(merged)
    cbc_to_seq = build_sequences(filtered)

    out_fa = Path(args.output_dir) / f"carrots_{args.sample}.fasta"
    n = write_fasta(cbc_to_seq, out_fa)
    print(f"[OK] {args.sample}: wrote {n} FASTA record(s) -> {out_fa}")


if __name__ == "__main__":
    main()

