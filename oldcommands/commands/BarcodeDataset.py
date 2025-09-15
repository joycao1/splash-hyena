import os, random, gzip, shutil, hashlib
from pathlib import Path
from typing import List, Optional, Union, Tuple, Dict
import torch
from pyfaidx import Fasta

class BarcodeGenomeDataset(torch.utils.data.Dataset):
    """
    Each FASTA record is a cell barcode 'genome'.
    Randomly sample a max_length window from a random barcode.
    Split is at the barcode level to avoid leakage.
    """
    def __init__(
        self,
        fasta_paths: List[str],          # one or more multi-record FASTAs
        split: str,                      # 'train' | 'valid' | 'test'
        max_length: int,
        total_size: int,
        tokenizer,
        tokenizer_name: str = "char",    # 'char' or 'bpe'
        add_eos: bool = False,
        pad_max_length: Optional[int] = None,  # used only for BPE
        task: str = "next_token_pred",   # or 'species_classification' -> class is file index
        auto_split_fracs: Tuple[float,float,float] = (0.8,0.1,0.1),
        auto_split_seed: int = 0,
        trim_fraction: float = 0.0,      # optionally trim ends of barcodes
        weights: str = "weighted_by_bp", # 'uniform' or 'weighted_by_bp'
    ):
        assert split in ("train","valid","test")
        self.split = split
        self.max_length = max_length
        self.total_size = total_size
        self.tokenizer = tokenizer
        self.tokenizer_name = tokenizer_name
        self.add_eos = add_eos
        self.pad_max_length = pad_max_length or max_length
        self.task = task

        # Open FASTAs (indexable)
        self.fastas: List[Fasta] = []
        self.records_per_file: List[List[str]] = []
        self.lengths_per_file: List[Dict[str,int]] = []

        for p in fasta_paths:
            p = Path(p)
            gz = p.with_suffix(p.suffix + ".gz")
            if gz.exists() and not p.exists():
                # One-time unzip if only .gz exists (optional—pyfaidx can also handle gz)
                with gzip.open(gz, "rb") as fin, open(p, "wb") as fout:
                    shutil.copyfileobj(fin, fout)

            fa = Fasta(str(p), sequence_always_upper=True)  # upper ⇒ tokenizer sees ACGTN
            self.fastas.append(fa)
            rec_ids = list(fa.keys())
            self.records_per_file.append(rec_ids)
            self.lengths_per_file.append({rid: len(fa[rid]) for rid in rec_ids})

        # Build deterministic barcode-level splits by hashing record ID
        tr, va, te = auto_split_fracs
        def bucket(bid: str) -> str:
            h = int(hashlib.blake2s((str(auto_split_seed)+bid).encode(), digest_size=4).hexdigest(), 16) / 2**32
            return "train" if h < tr else ("valid" if h < tr+va else "test")

        self.pool: List[Tuple[int,str,int]] = []  # (file_idx, record_id, length)
        for fi, rec_ids in enumerate(self.records_per_file):
            for rid in rec_ids:
                if bucket(rid) == split:
                    self.pool.append((fi, rid, self.lengths_per_file[fi][rid]))

        # Sampling weights across barcodes
        if weights == "uniform":
            self.pool_weights = None
        elif weights == "weighted_by_bp":
            tot = sum(L for _,_,L in self.pool) or 1
            self.pool_weights = [L/tot for _,_,L in self.pool]
        else:
            raise ValueError("weights must be 'uniform' or 'weighted_by_bp'")

        # For classification, label by file index (or customize)
        self.d_output = len(self.fastas) if self.task == "species_classification" else None
        self.trim_fraction = trim_fraction

    def __len__(self):
        return self.total_size

    def __getitem__(self, idx):
        fi, rid, L = random.choices(self.pool, weights=self.pool_weights, k=1)[0]
        rec = self.fastas[fi][rid]

        # Allowed start range (optionally trim ends), and ensure a full window fits
        trim = int(L * self.trim_fraction)
        max_start = max(0, (L - self.max_length) - trim)
        left = trim
        right_start = max(left, max_start)
        start = random.randint(left, right_start) if right_start >= left else 0
        end = start + self.max_length

        seq = str(rec[start:end])
        if len(seq) < self.max_length:
            # pad RIGHT with 'N' (semantic unknown base; keep model pad '4' in token space if needed)
            seq = seq.ljust(self.max_length, "N")

        # Tokenize
        if self.tokenizer_name == "char":
            tok = self.tokenizer(seq, add_special_tokens=False)["input_ids"]
            if self.add_eos:
                tok.append(self.tokenizer.sep_token_id)
        elif self.tokenizer_name == "bpe":
            out = self.tokenizer(seq, padding="max_length", max_length=self.pad_max_length, truncation=True)
            ids = out["input_ids"]
            tok = ids[1:] if self.add_eos else ids[1:-1]  # mirror your earlier logic
        else:
            raise ValueError("tokenizer_name must be 'char' or 'bpe'")

        seq_ids = torch.LongTensor(tok)
        data = seq_ids[:-1].clone()
        if self.task == "next_token_pred":
            target = seq_ids[1:].clone()
        elif self.task == "species_classification":
            target = fi
        else:
            raise ValueError("task must be 'next_token_pred' or 'species_classification'")

        return data, target

