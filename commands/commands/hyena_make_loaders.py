#!/usr/bin/env python3
import argparse, glob, os
from typing import List
import torch
from torch.utils.data import DataLoader

# Your minimal barcode-as-genome dataset
from bar_gene_dataset import BarcodeGenomeDataset

# Your tokenizer (char-level for A/C/G/T/N/'4')
from src.dataloaders.datasets.hg38_char_tokenizer import CharacterTokenizer

def _paths_from_glob(pattern: str) -> List[str]:
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No files matched: {pattern}")
    return paths

def build_loaders(
    train_glob: str,
    valid_glob: str = "",
    test_glob: str = "",
    max_length: int = 250_000,
    total_size_train: int = 200_000,
    total_size_valid: int = 10_000,
    total_size_test: int = 10_000,
    batch_size: int = 1,
    num_workers_train: int = 4,
    num_workers_eval: int = 2,
    tokenizer_lowercase: bool = False,
):
    # Tokenizer
    tokenizer = CharacterTokenizer()

    # Resolve shards
    train_paths = _paths_from_glob(train_glob)
    valid_paths = _paths_from_glob(valid_glob) if valid_glob else train_paths
    test_paths  = _paths_from_glob(test_glob)  if test_glob  else train_paths

    # If you only have one pool of FASTAs, we’ll auto-split by barcode ID.
    def make_ds(paths, split, total_size):
        # If you pre-split by different shard sets per split, set auto_split_fracs=(1,0,0)
        auto = (0.8, 0.1, 0.1) if (train_paths == valid_paths == test_paths) else \
               (1.0, 0.0, 0.0) if split == "train" else \
               (0.0, 1.0, 0.0) if split == "valid" else (0.0, 0.0, 1.0)

        return BarcodeGenomeDataset(
            fasta_paths=paths,
            split=split,
            max_length=max_length,
            total_size=total_size,
            tokenizer=tokenizer,
            tokenizer_name="char",
            add_eos=False,
            task="next_token_pred",
            auto_split_fracs=auto,
            auto_split_seed=0,
            trim_fraction=0.0,
            weights="weighted_by_bp",
        )

    train_ds = make_ds(train_paths, "train", total_size_train)
    valid_ds = make_ds(valid_paths, "valid", total_size_valid)
    test_ds  = make_ds(test_paths,  "test",  total_size_test)

    # DataLoaders (don’t shuffle; dataset already samples randomly)
    train_loader = DataLoader(
        train_ds, batch_size=batch_size, shuffle=False,
        num_workers=num_workers_train, prefetch_factor=2,
        persistent_workers=(num_workers_train > 0), pin_memory=True
    )
    valid_loader = DataLoader(
        valid_ds, batch_size=batch_size, shuffle=False,
        num_workers=num_workers_eval, prefetch_factor=2,
        persistent_workers=(num_workers_eval > 0), pin_memory=True
    )
    test_loader = DataLoader(
        test_ds, batch_size=batch_size, shuffle=False,
        num_workers=num_workers_eval, prefetch_factor=2,
        persistent_workers=(num_workers_eval > 0), pin_memory=True
    )

    return train_loader, valid_loader, test_loader

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_glob", required=True, help="e.g. fasta_data/*.fasta")
    ap.add_argument("--valid_glob", default="", help="optional; else auto-split")
    ap.add_argument("--test_glob",  default="", help="optional; else auto-split")
    ap.add_argument("--max_length", type=int, default=250_000)
    ap.add_argument("--total_size_train", type=int, default=200_000)
    ap.add_argument("--total_size_valid", type=int, default=10_000)
    ap.add_argument("--total_size_test",  type=int, default=10_000)
    ap.add_argument("--batch_size", type=int, default=1)
    ap.add_argument("--num_workers_train", type=int, default=4)
    ap.add_argument("--num_workers_eval", type=int, default=2)
    args = ap.parse_args()

    tl, vl, _ = build_loaders(
        train_glob=args.train_glob,
        valid_glob=args.valid_glob,
        test_glob=args.test_glob,
        max_length=args.max_length,
        total_size_train=args.total_size_train,
        total_size_valid=args.total_size_valid,
        total_size_test=args.total_size_test,
        batch_size=args.batch_size,
        num_workers_train=args.num_workers_train,
        num_workers_eval=args.num_workers_eval,
    )

    # Smoke test
    x, y = next(iter(tl))
    print("Train batch shapes:", tuple(x.shape), tuple(y.shape))
    vx, vy = next(iter(vl))
    print("Valid batch shapes:", tuple(vx.shape), tuple(vy.shape))
    # NOTE: You normally import build_loaders() from your train script.

if __name__ == "__main__":
    main()

