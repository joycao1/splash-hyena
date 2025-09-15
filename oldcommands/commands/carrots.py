#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert BKC-style TSV -> concatenated-per-CBC FASTA, after removing targets
that appear across all anchors (evaluated per CBC).

TSV columns (header or --no-header):
  sample_id, cbc, anchor, target, count

Process:
  1) Sort by (cbc, anchor, count desc)
  2) For each CBC, drop any target that appears in *every* anchor (only when CBC has >1 anchors)
  3) For each CBC, concatenate: [anchor] + all [targets for that anchor in sorted order]
  4) Emit one FASTA record per CBC

Example:
  python tsv_to_fasta_filtered.py -i data.tsv -o ./fasta_out --output-name carrots.fasta
"""

import argparse
from pathlib import Path
import sys

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def parse_args():
    ap = argparse.ArgumentParser(
        description="TSV -> FASTA with per-CBC universal-target filtering."
    )
    ap.add_argument("-i", "--input", required=True, help="Path to input TSV file.")
    ap.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory to write outputs (created if missing).",
    )
    ap.add_argument(
        "--output-name",
        default="output.fasta",
        help="Output FASTA filename (default: output.fasta).",
    )
    ap.add_argument(
        "--no-header",
        action="store_true",
        help="Set if the TSV has no header row.",
    )
    ap.add_argument(
        "--print-head",
        type=int,
        default=5,
        help="Rows to print from head for a quick check (0 disables).",
    )
    return ap.parse_args()


def sanitize_fasta_id(s: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(ch if ch in allowed else "_" for ch in str(s))


def load_tsv(path: Path, no_header: bool) -> pd.DataFrame:
    names = ["sample_id", "cbc", "anchor", "target", "count"]
    read_kwargs = dict(sep="\t")
    if no_header:
        read_kwargs.update(header=None, names=names)

    df = pd.read_csv(path, **read_kwargs)

    # Validate required cols
    required = ["cbc", "anchor", "target", "count"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")

    # Normalize types
    df["cbc"] = df["cbc"].astype(str)
    df["anchor"] = df["anchor"].astype(str)
    df["target"] = df["target"].astype(str)
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)

    # Drop empties
    df = df[(df["anchor"].str.len() > 0) & (df["target"].str.len() > 0)].copy()

    return df


def filter_universal_targets(df_sorted: pd.DataFrame) -> pd.DataFrame:
    """
    Drop rows where (cbc, target) appears in *all* anchors for that CBC.
    Only applies when CBC has >1 distinct anchors; otherwise keep everything.
    """
    # number of distinct anchors per CBC
    n_anchors_per_cbc = df_sorted.groupby("cbc")["anchor"].nunique().rename("n_anchors")

    # for each (cbc, target), how many distinct anchors does it show up in?
    anchors_per_target = (
        df_sorted.groupby(["cbc", "target"])["anchor"]
        .nunique()
        .rename("anchors_for_target")
        .reset_index()
    )

    df_tmp = df_sorted.merge(anchors_per_target, on=["cbc", "target"], how="left")
    df_tmp = df_tmp.merge(n_anchors_per_cbc, on="cbc", how="left")

    # keep if NOT universal, or if only one anchor exists
    keep_mask = ~(
        (df_tmp["anchors_for_target"] == df_tmp["n_anchors"])
        & (df_tmp["n_anchors"] > 1)
    )

    df_filtered = df_tmp.loc[keep_mask].drop(columns=["anchors_for_target", "n_anchors"])
    return df_filtered


def build_sequences(df_filtered: pd.DataFrame) -> dict:
    """
    Build concatenated sequence per CBC:
      for each CBC:
        for each anchor (post-sort order):
          append anchor once, then all targets for that anchor (post-sort order)
    """
    cbc_to_seq = {}
    for cbc, group in df_filtered.groupby("cbc", sort=False):
        parts = []
        for anchor, a_group in group.groupby("anchor", sort=False):
            parts.append(anchor)
            parts.extend(a_group["target"].tolist())
        cbc_to_seq[cbc] = "".join(parts)
    return cbc_to_seq


def write_fasta(cbc_to_seq: dict, out_path: Path) -> int:
    records = [
        SeqRecord(Seq(seq), id=sanitize_fasta_id(cbc), description="")
        for cbc, seq in cbc_to_seq.items()
    ]
    return SeqIO.write(records, str(out_path), "fasta")


def main():
    args = parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        print(f"ERROR: Input file not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / args.output_name

    df = load_tsv(in_path, args.no_header)

    if args.print_head and args.print_head > 0:
        print(df.head(args.print_head).to_string(index=False))

    # Sort by count (descending) for each CBC/anchor pair
    df_sorted = df.sort_values(
        by=["cbc", "anchor", "count"], ascending=[True, True, False]
    ).copy()

    # Drop per-CBC universal targets
    df_filtered = filter_universal_targets(df_sorted)

    # Build concatenated sequences
    cbc_to_seq = build_sequences(df_filtered)

    # Write FASTA
    n = write_fasta(cbc_to_seq, out_path)
    print(f"Wrote {n} FASTA record(s) to: {out_path}")


if __name__ == "__main__":
    main()
