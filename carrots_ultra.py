#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ultra-low-memory per-sample merge -> FASTA.

Pipeline (outside this script):
  cat chunks -> normalize to 4 cols -> sort+sum on disk
This script:
  - reads aggregated lines sorted by cbc,anchor,count(desc),
  - filters universal targets per CBC,
  - writes FASTA incrementally.

Input format (one per line, tab-separated):
  cbc <TAB> anchor <TAB> target <TAB> count
MUST be sorted by: cbc (asc), anchor (asc), count (desc within anchor).
"""

import sys
from pathlib import Path

def sanitize_id(s: str) -> str:
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-")
    return "".join(ch if ch in allowed else "_" for ch in s)

def write_fasta_record(w, cbc, anchors_targets):
    """anchors_targets: list of (anchor, [targets_in_order]) for this CBC"""
    parts = []
    for anchor, targets in anchors_targets:
        parts.append(anchor)
        parts.extend(targets)
    seq = "".join(parts)
    w.write(f">{sanitize_id(cbc)}\n{seq}\n")

def process_sorted_agg(agg_path: Path, out_fasta: Path):
    with agg_path.open("r") as f, out_fasta.open("w") as w:
        cur_cbc = None
        # For universal filtering within a CBC:
        anchors_for_cbc = set()
        targets_to_anchors = dict()  # target -> set(anchors)
        anchor_to_targets = []       # list of (anchor, [targets]) preserving order

        def flush():
            nonlocal cur_cbc, anchors_for_cbc, targets_to_anchors, anchor_to_targets
            if cur_cbc is None:
                return
            n_anchors = len(anchors_for_cbc) if anchors_for_cbc else 0
            # Filter out targets that appear in ALL anchors for this CBC (only if >1 anchors)
            filtered = []
            for anchor, targets in anchor_to_targets:
                kept = []
                for t in targets:
                    s = targets_to_anchors.get(t, set())
                    if not (n_anchors > 1 and len(s) == n_anchors):
                        kept.append(t)
                filtered.append((anchor, kept))
            # Write FASTA record
            write_fasta_record(w, cur_cbc, filtered)
            # reset
            cur_cbc = None
            anchors_for_cbc = set()
            targets_to_anchors = dict()
            anchor_to_targets = []

        # Read aggregated, sorted lines
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 4:
                continue
            cbc, anchor, target, cnt = parts
            # cnt is not used now (already sorted by count), presence is enough
            if cbc != cur_cbc:
                flush()
                cur_cbc = cbc
                anchors_for_cbc = set()
                targets_to_anchors = {}
                anchor_to_targets = []

            if not anchor_to_targets or anchor_to_targets[-1][0] != anchor:
                anchor_to_targets.append((anchor, []))
            anchor_to_targets[-1][1].append(target)
            anchors_for_cbc.add(anchor)
            s = targets_to_anchors.get(target)
            if s is None:
                s = set()
                targets_to_anchors[target] = s
            s.add(anchor)

        flush()

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} AGG_SORTED_TSV OUT_FASTA", file=sys.stderr)
        sys.exit(2)
    agg = Path(sys.argv[1])
    out_fa = Path(sys.argv[2])
    process_sorted_agg(agg, out_fa)
    print(f"[OK] wrote FASTA -> {out_fa}", file=sys.stderr)

if __name__ == "__main__":
    main()

