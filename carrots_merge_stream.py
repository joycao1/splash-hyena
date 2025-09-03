#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Streaming merge of bkctxt chunks into one FASTA per sample, low memory.

Input files (per sample): bkctxt/<SAMPLE>_R1.part_XX.txt
Each line has 4 or 5 whitespace-separated columns:
  4-col: cbc   anchor   target   count
  5-col: <ignored> cbc  anchor   target   count

We:
  1) Sum counts across all chunks for identical (cbc, anchor, target).
  2) Drop targets that appear in ALL anchors for a CBC (only if CBC has >1 anchors).
  3) For each CBC, concatenate: [anchor] + all [targets] for that anchor, with targets sorted by count desc.
  4) Write one FASTA per sample to out/fasta/carrots_<SAMPLE>.fasta
"""

from pathlib import Path
from collections import defaultdict, Counter
import argparse
import sys

def parse_args():
    ap = argparse.ArgumentParser(description="Streaming per-sample merge to FASTA (low memory).")
    ap.add_argument("-s", "--sample", required=True, help="Sample ID, e.g. ERR13720419")
    ap.add_argument("-i", "--input-dir", default="bkctxt", help="Dir with chunk .txt files")
    ap.add_argument("-o", "--output-dir", default="out/fasta", help="Dir to write FASTA")
    ap.add_argument("--print-head", type=int, default=0, help="Print first N merged rows (stderr)")
    return ap.parse_args()

def sanitize_id(s: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(ch if ch in allowed else "_" for ch in s)

def stream_sum_counts(sample: str, input_dir: Path):
    """
    Return:
      counts: dict[(cbc, anchor, target)] -> int
      order_hint: dict[cbc] -> dict[anchor] -> total_count (used for stable-ish ordering)
    """
    counts = defaultdict(int)
    order_hint = defaultdict(lambda: defaultdict(int))
    files = sorted((input_dir).glob(f"{sample}_R1.part_*.txt"))
    if not files:
        raise FileNotFoundError(f"No input files matched {input_dir}/{sample}_R1.part_*.txt")
    for fp in files:
        with fp.open("r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()  # whitespace
                if len(parts) == 5:
                    _, cbc, anchor, target, cnt = parts
                elif len(parts) == 4:
                    cbc, anchor, target, cnt = parts
                else:
                    # skip malformed lines gently
                    continue
                try:
                    c = int(cnt)
                except ValueError:
                    continue
                key = (cbc, anchor, target)
                counts[key] += c
                order_hint[cbc][anchor] += c
    return counts, order_hint

def filter_universal(counts):
    """
    Remove (cbc,target) that appear in ALL anchors for that cbc (if >1 anchors).
    Return new dict with same shape.
    """
    # map: cbc -> set(anchors), and cbc->target->anchors_count
    anchors_by_cbc = defaultdict(set)
    target_anchor_count = defaultdict(lambda: defaultdict(int))  # cbc -> target -> n_anchors
    presence = defaultdict(lambda: defaultdict(set))  # cbc -> target -> set(anchors)

    for (cbc, anchor, target), c in counts.items():
        anchors_by_cbc[cbc].add(anchor)
        presence[cbc][target].add(anchor)

    for cbc, targ_map in presence.items():
        for target, anchors in targ_map.items():
            target_anchor_count[cbc][target] = len(anchors)

    # build filtered
    filtered = {}
    for (cbc, anchor, target), c in counts.items():
        n_anchors = len(anchors_by_cbc[cbc])
        if n_anchors > 1 and target_anchor_count[cbc][target] == n_anchors:
            # universal across anchors for this CBC -> drop
            continue
        filtered[(cbc, anchor, target)] = c
    return filtered

def build_sequences(counts, order_hint):
    """
    Build concatenated sequence per CBC.
    Anchor order: by descending total anchor count (from order_hint).
    Targets within anchor: by descending count.
    """
    by_cbc_anchor = defaultdict(lambda: defaultdict(list))  # cbc -> anchor -> [(target, count)]
    for (cbc, anchor, target), c in counts.items():
        by_cbc_anchor[cbc][anchor].append((target, c))

    cbc_to_seq = {}
    for cbc, anchor_map in by_cbc_anchor.items():
        # order anchors by total count desc; fallback to name
        anchors_sorted = sorted(
            anchor_map.keys(),
            key=lambda a: (-order_hint[cbc].get(a, 0), a)
        )
        parts = []
        for anchor in anchors_sorted:
            parts.append(anchor)
            # sort targets by count desc, then lexicographically for stability
            targets_sorted = sorted(anchor_map[anchor], key=lambda tc: (-tc[1], tc[0]))
            parts.extend([t for (t, _) in targets_sorted])
        cbc_to_seq[cbc] = "".join(parts)
    return cbc_to_seq

def write_fasta(cbc_to_seq, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as w:
        for cbc, seq in cbc_to_seq.items():
            w.write(f">{sanitize_id(cbc)}\n{seq}\n")

def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    counts, order_hint = stream_sum_counts(args.sample, input_dir)
    # optional preview: print some merged rows
    if args.print_head > 0:
        # build a small preview without materializing everything
        shown = 0
        for (cbc, anchor, target), c in counts.items():
            print(f"{cbc}\t{anchor}\t{target}\t{c}", file=sys.stderr)
            shown += 1
            if shown >= args.print_head:
                break

    filtered = filter_universal(counts)
    cbc_to_seq = build_sequences(filtered, order_hint)
    out_fa = output_dir / f"carrots_{args.sample}.fasta"
    write_fasta(cbc_to_seq, out_fa)
    print(f"[OK] {args.sample}: wrote {len(cbc_to_seq)} FASTA record(s) -> {out_fa}")

if __name__ == "__main__":
    main()

