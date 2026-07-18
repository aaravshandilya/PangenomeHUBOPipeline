#!/usr/bin/env python3
"""
Stage 1c - Problem creation: simulate shotgun short-read sequencing.

Single-ended reads, random strand, uniform substitution errors, no indel
errors (matches Methods -> "Annotating the graph with copy numbers").

Usage:
    shred_genome.py --fasta test_genomes.fasta --config config.yaml \
                     --seed 1 --out-dir results/seed_1/reads
"""
import argparse
import random
from pathlib import Path

import yaml

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def read_fasta(path):
    records = []
    name, seq = None, []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def simulate_reads(rng, genome, coverage, read_len, err_rate):
    n_reads = max(1, (len(genome) * coverage) // read_len)
    reads = []
    for i in range(n_reads):
        if len(genome) <= read_len:
            start = 0
        else:
            start = rng.randint(0, len(genome) - read_len)
        frag = genome[start:start + read_len]
        if rng.random() < 0.5:
            frag = revcomp(frag)
        frag = list(frag)
        for j in range(len(frag)):
            if rng.random() < err_rate:
                frag[j] = rng.choice([b for b in BASES if b != frag[j]])
        reads.append((f"read_{i}_pos{start}", "".join(frag)))
    return reads


def write_fastq(path, reads):
    with open(path, "w") as fh:
        for name, seq in reads:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))["sequencing"]
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed + 9999)
    for name, seq in read_fasta(args.fasta):
        reads = simulate_reads(rng, seq, cfg["coverage"], cfg["read_length"],
                                cfg["substitution_error_rate"])
        write_fastq(out / f"{name}.fastq", reads)
        print(f"[shred_genome] {name}: {len(reads)} reads "
              f"({cfg['coverage']}x coverage, {cfg['read_length']}bp)")


if __name__ == "__main__":
    main()
