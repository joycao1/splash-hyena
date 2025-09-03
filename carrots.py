#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
carrots.py  —  Drop-in replacement

Two modes:

A) Single-file mode (backwards-compatible)
   -i INPUT [-o OUTDIR] [--output-name NAME] [--no-header] [--print-head N]
   INPUT is a 4- or 5-column TSV/whitespace file:
      4-col: cbc  anchor  target  count
      5-col: IGNORE  cbc  anchor  target  count   (first column often '0')
   We group by (cbc, anchor, target), SUM counts (handles duplicates within the file),
   drop "universal" targets (target present in *all* anchors for a CBC when CBC has >1 anchors),
   then build one FASTA record per CBC by concatenating [anchor] + all [targets].

B) Per-sample merge mode
   --sample SAMPLE [-i INPUT_DIR] [-o OUTDIR] [--output-name NAME] [--print-head N]
   Reads many chunk files: INPUT_DIR/<SAMPLE>_R1.part_*.txt (default INPUT_DIR='bkctxt'),
   sums counts across chunks, applies the same filter, writes one FASTA:
      OUTDIR/(output-name or 'carrots_<SAMPLE>.fasta')

No heavy deps; streaming I/O; low memory.
"""

from pathlib import Path
from collections import defaultdict
import argparse
import sys
import re
from typing import Iterable, Tuple, Dict, List, Optional

# ---------- CLI ----------

def parse_args():
    ap = argparse.ArgumentParser(description="TSV/chunks -> FASTA with universal-target filtering (low memory).")
    m = ap.add_mutually_exclusive_group(required=True)
    m.add_argument("-i", "--input",
                   help="Single input file (4 or 5 cols). In sample mode, this is INPUT_DIR (default: bkctxt).")
    m.add_argument("--sample",
                   help="Sample ID for per-sample merge mode, matches <S>_R1.part_*.txt in INPUT_DIR (default: bkctxt).")

    ap.add_argument("-o", "--output-dir", default="out/fasta", help="Directory for FASTA output.")
    ap.add_argument("--output-name", default=None,
                    help="Output FASTA filename. If unset: single-file mode -> 'output.fasta'; sample mode -> 'carrots_<S>.fasta'.")
    ap.add_argument("--no-header", action="store_true",
                    help="Set only if single-file INPUT has no header row. (Sample mode ignores header anyway.)")
    ap.add_argument("--print-head", type=int, default=0,
                    help="Print first N merged rows (cbc\\tanchor\\ttarget\\tcount) to stderr for sanity.")
    return ap.parse_args()

# ---------- Parsing & merging (streaming) ----------

def _split_line(line: str) -> List[str]:
    """Split on tabs if present, else on any whitespace."""
    line = line.rstrip("\n\r")
    if "\t" in line:
        return line.split("\t")
    return line.split()

def _read_rows_stream(fp: Iterable[str],
                      expect_header: Optional[bool]) -> Iterable[Tuple[str, str, str, int]]:
    """
    Yield (cbc, anchor, target, count) from a file that has 4 or 5 columns.
    - If expect_header is True, drop the first line.
    - If expect_header is False or None, try to auto-skip a header if the last column isn't an int.
    """
    first = True
    for raw in fp:
        raw = raw.strip()
        if not raw:
            continue
        parts = _split_line(raw)
        # 5 or 4 columns only
        if len(parts) == 5:
            parts = parts[1:]  # drop first
        elif len(parts) != 4:
            # ignore malformed lines gracefully
            continue

        cbc, anchor, target, cnt = parts

        # header handling on first data row
        if first:
            first = False
            if expect_header is True:
                # drop this line as header and continue
                # (but only if it looks non-numeric in count column)
                try:
                    int(cnt)
                except ValueError:
                    # skip header
                    continue
            elif expect_header is None:
                # auto-detect: if count not int, skip as header
                try:
                    int(cnt)
                except ValueError:
                    continue

        # convert count
        try:
            c = int(cnt)
        except ValueError:
            # non-numeric count → skip
            continue

        if not anchor or not target:
            continue

        yield (str(cbc), str(anchor), str(target), c)

def merge_file(path: Path, no_header: bool) -> Dict[Tuple[str, str, str], int]:
    """Read a single file and sum counts for duplicate (cbc, anchor, target) rows."""
    counts = defaultdict(int)
    expect_header = (False if no_header else None)  # None = auto-detect
    with path.open("r") as f:
        for cbc, anchor, target, c in _read_rows_stream(f, expect_header):
            counts[(cbc, anchor, target)] += c
    return counts

def merge_sample_chunks(sample: str, input_dir: Path) -> Dict[Tuple[str, str, str], int]:
    """Read all chunks for a sample and sum counts across files."""
    counts = defaultdict(int)
    pattern = f"{sample}_R1.part_*.txt"
    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matched: {input_dir}/{pattern}")
    for path in files:
        with path.open("r") as f:
            for cbc, anchor, target, c in _read_rows_stream(f, expect_header=False):
                counts[(cbc, anchor, target)] += c
    return counts

# ---------- Filtering & sequence building ----------

def filter_universal(counts: Dict[Tuple[str, str, str], int]) -> Dict[Tuple[str, str, str], int]:
    """
    Remove rows where (cbc, target) appears in ALL anchors for that cbc
    (only when CBC has >1 anchors).
    """
    anchors_by_cbc = defaultdict(set)
    targ_anchors = defaultdict(lambda: defaultdict(set))  # cbc -> target -> set(anchors)

    for (cbc, anchor, target), _c in counts.items():
        anchors_by_cbc[cbc].add(anchor)
        targ_anchors[cbc][target].add(anchor)

    filtered = {}
    for (cbc, anchor, target), c in counts.items():
        n_anchors = len(anchors_by_cbc[cbc])
        if n_anchors > 1 and len(targ_anchors[cbc][target]) == n_anchors:
            continue
        filtered[(cbc, anchor, target)] = c
    return filtered

def build_sequences(counts: Dict[Tuple[str, str, str], int]) -> Dict[str, str]:
    """
    For each CBC:
      - order anchors by total count desc (tie-break lexicographically),
      - within each anchor, order targets by count desc then lexicographically,
      - concatenate anchor + targets.
    """
    anchor_totals = defaultdict(lambda: defaultdict(int))   # cbc -> anchor -> total count
    by_cbc_anchor = defaultdict(lambda: defaultdict(list))  # cbc -> anchor -> [(target, count)]

    for (cbc, anchor, target), c in counts.items():
        anchor_totals[cbc][anchor] += c
        by_cbc_anchor[cbc][anchor].append((target, c))

    cbc_to_seq = {}
    for cbc, amap in by_cbc_anchor.items():
        anchors_sorted = sorted(amap.keys(), key=lambda a: (-anchor_totals[cbc][a], a))
        parts: List[str] = []
        for anchor in anchors_sorted:
            parts.append(anchor)
            t_sorted = sorted(amap[anchor], key=lambda tc: (-tc[1], tc[0]))
            parts.extend([t for (t, _) in t_sorted])
        cbc_to_seq[cbc] = "".join(parts)
    return cbc_to_seq

# ---------- FASTA output & utils ----------

def _sanitize_id(s: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(ch if ch in allowed else "_" for ch in s)

def write_fasta(cbc_to_seq: Dict[str, str], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as w:
        for cbc, seq in cbc_to_seq.items():
            w.write(f">{_sanitize_id(cbc)}\n{seq}\n")

def print_head_preview(counts: Dict[Tuple[str, str, str], int], n: int) -> None:
    shown = 0
    for (cbc, anchor, target), c in counts.items():
        print(f"{cbc}\t{anchor}\t{target}\t{c}", file=sys.stderr)
        shown += 1
        if shown >= n:
            break

# ---------- Main ----------

def main():
    args = parse_args()
    out_dir = Path(args.output_dir)

    if args.sample:
        # per-sample merge mode
        sample = args.sample
        input_dir = Path(args.input or "bkctxt")
        counts = merge_sample_chunks(sample, input_dir)
        if args.print_head > 0:
            print_head_preview(counts, args.print_head)
        counts = filter_universal(counts)
        cbc_to_seq = build_sequences(counts)
        out_name = args.output_name or f"carrots_{sample}.fasta"
        out_path = out_dir / out_name
        write_fasta(cbc_to_seq, out_path)
        print(f"[OK] {sample}: wrote {len(cbc_to_seq)} FASTA record(s) -> {out_path}")
    else:
        # single-file mode (back-compat)
        in_path = Path(args.input)
        if not in_path.exists():
            print(f"ERROR: Input file not found: {in_path}", file=sys.stderr)
            sys.exit(1)
        counts = merge_file(in_path, no_header=args.no_header)
        if args.print_head > 0:
            print_head_preview(counts, args.print_head)
        counts = filter_universal(counts)
        cbc_to_seq = build_sequences(counts)
        out_name = args.output_name or "output.fasta"
        out_path = out_dir / out_name
        write_fasta(cbc_to_seq, out_path)
        print(f"[OK] wrote {len(cbc_to_seq)} FASTA record(s) -> {out_path}")

if __name__ == "__main__":
    main()

