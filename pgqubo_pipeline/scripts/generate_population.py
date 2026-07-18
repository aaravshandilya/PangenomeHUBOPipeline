#!/usr/bin/env python3
"""
Stage 1a - Problem creation: synthetic genome population.

Mirrors the paper's "genome_create" tool (Methods -> Pangenome creation):
    - a founder genome with STRs, CNVs, large repeats, translocations,
      inversions and point mutations
    - 10 generations of haploid descent, each derived from a single
      random parent and lightly re-mutated
    - a random subset used to build the pangenome, a disjoint subset
      held out as "new individuals" to be assembled

Usage:
    generate_population.py --config config.yaml --seed 1 --out-dir results/seed_1/population
"""
import argparse
import json
import random
from pathlib import Path

import numpy as np
import yaml

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def random_seq(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def insert_str(rng, seq, unit_len_choices=(2, 3, 4, 5)):
    unit = random_seq(rng, rng.choice(unit_len_choices))
    repeats = rng.randint(5, 20)
    pos = rng.randint(0, len(seq))
    return seq[:pos] + unit * repeats + seq[pos:]


def apply_cnv(rng, seq, min_len=200, max_len=2000):
    """Duplicate or delete a chunk of sequence."""
    length = rng.randint(min_len, max_len)
    if len(seq) <= length + 1:
        return seq
    pos = rng.randint(0, len(seq) - length)
    chunk = seq[pos:pos + length]
    if rng.random() < 0.5:  # duplication, inserted nearby
        insert_at = rng.randint(0, len(seq))
        return seq[:insert_at] + chunk + seq[insert_at:]
    else:  # deletion
        return seq[:pos] + seq[pos + length:]


def apply_large_repeat(rng, seq, min_len=500, max_len=3000):
    length = rng.randint(min_len, max_len)
    if len(seq) <= length + 1:
        return seq
    pos = rng.randint(0, len(seq) - length)
    chunk = seq[pos:pos + length]
    insert_at = rng.randint(0, len(seq))
    return seq[:insert_at] + chunk + seq[insert_at:]


def apply_translocation(rng, seq, min_len=200, max_len=1500):
    length = rng.randint(min_len, max_len)
    if len(seq) <= 2 * length + 2:
        return seq
    pos = rng.randint(0, len(seq) - length)
    chunk = seq[pos:pos + length]
    remainder = seq[:pos] + seq[pos + length:]
    insert_at = rng.randint(0, len(remainder))
    return remainder[:insert_at] + chunk + remainder[insert_at:]


def apply_inversion(rng, seq, min_len=200, max_len=1500):
    length = rng.randint(min_len, max_len)
    if len(seq) <= length + 1:
        return seq
    pos = rng.randint(0, len(seq) - length)
    chunk = seq[pos:pos + length]
    return seq[:pos] + revcomp(chunk) + seq[pos + length:]


def apply_point_mutations(rng, seq, rate):
    seq = list(seq)
    n_mut = np.random.binomial(len(seq), rate)
    for _ in range(n_mut):
        i = rng.randrange(len(seq))
        seq[i] = rng.choice([b for b in BASES if b != seq[i]])
    return "".join(seq)


def mutate_genome(rng, seq, rates, multiplier=1.0):
    if rng.random() < rates["str_indel"] * multiplier:
        seq = insert_str(rng, seq)
    if rng.random() < rates["cnv"] * multiplier:
        seq = apply_cnv(rng, seq)
    if rng.random() < rates["large_repeat"] * multiplier:
        seq = apply_large_repeat(rng, seq)
    if rng.random() < rates["translocation"] * multiplier:
        seq = apply_translocation(rng, seq)
    if rng.random() < rates["inversion"] * multiplier:
        seq = apply_inversion(rng, seq)
    seq = apply_point_mutations(rng, seq, rates["point"] * multiplier)
    return seq


def build_population(cfg, seed):
    rng = random.Random(seed)
    np.random.seed(seed)

    pop_cfg = cfg["population"]
    rates = pop_cfg["mutation_rates"]

    founder = random_seq(rng, pop_cfg["genome_length"])
    founder = mutate_genome(rng, founder, rates, multiplier=pop_cfg["root_mutation_multiplier"])

    genomes = {0: founder}
    lineage = {0: None}

    n_members = pop_cfg["n_members"]
    n_generations = pop_cfg["n_generations"]
    per_gen = max(1, (n_members - 1) // n_generations)

    next_id = 1
    current_gen_ids = [0]
    for gen in range(n_generations):
        new_gen_ids = []
        n_this_gen = per_gen if gen < n_generations - 1 else (n_members - next_id)
        for _ in range(max(0, n_this_gen)):
            if next_id >= n_members:
                break
            parent_id = rng.choice(current_gen_ids)
            child_seq = mutate_genome(rng, genomes[parent_id], rates, multiplier=1.0)
            genomes[next_id] = child_seq
            lineage[next_id] = parent_id
            new_gen_ids.append(next_id)
            next_id += 1
        current_gen_ids = new_gen_ids or current_gen_ids

    all_ids = list(genomes.keys())
    rng.shuffle(all_ids)
    pangenome_ids = sorted(all_ids[:pop_cfg["n_pangenome_build"]])
    remaining = [i for i in all_ids if i not in pangenome_ids]
    test_ids = sorted(remaining[:cfg["n_test_genomes"]])

    return genomes, lineage, pangenome_ids, test_ids


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    genomes, lineage, pangenome_ids, test_ids = build_population(cfg, args.seed)

    write_fasta(out / "pangenome_build_genomes.fasta",
                [(f"member_{i}", genomes[i]) for i in pangenome_ids])
    write_fasta(out / "test_genomes.fasta",
                [(f"test_{i}", genomes[i]) for i in test_ids])

    meta = {
        "seed": args.seed,
        "n_members": len(genomes),
        "pangenome_ids": pangenome_ids,
        "test_ids": test_ids,
        "lineage": lineage,
    }
    json.dump(meta, open(out / "population_meta.json", "w"), indent=2)
    print(f"[generate_population] seed={args.seed}: "
          f"{len(genomes)} genomes, {len(pangenome_ids)} for pangenome, "
          f"{len(test_ids)} held out for testing")


if __name__ == "__main__":
    main()
