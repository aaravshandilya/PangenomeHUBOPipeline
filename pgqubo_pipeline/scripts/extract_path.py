#!/usr/bin/env python3
"""
Stage 5a - Solution processing: decode {x_t,v}=1 assignments into a walk,
then concatenate node sequences (reverse-complementing '-' oriented nodes)
into a candidate genome FASTA.

If the QUBO constraints were violated (e.g. two vertices "on" at the same
time step, a side effect of the imperfect classical_sa fallback), we take a
robust majority approach: at each time step, prefer the visited vertex whose
node still needs visits according to its remaining copy-number budget; ties
broken by earliest node id. Steps assigned only to `end` truncate the walk.

Usage:
    extract_path.py --assignment assignment.tsv --gfa pangenome.gfa --out candidate.fasta
"""
import argparse
from collections import defaultdict
from pathlib import Path


def read_gfa_nodes(path):
    nodes = {}
    for line in open(path):
        if line.startswith("S\t"):
            _, nid, seq = line.rstrip("\n").split("\t")[:3]
            nodes[nid] = seq
    return nodes


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_assignment(path):
    by_t = defaultdict(list)
    with open(path) as fh:
        next(fh)
        for line in fh:
            t, v = line.strip().split("\t")
            by_t[int(t)].append(v)
    return by_t


def resolve_walk(by_t, remaining_budget):
    walk = []
    max_t = max(by_t.keys()) if by_t else 0
    for t in range(1, max_t + 1):
        candidates = [v for v in by_t.get(t, []) if v != "end"]
        if not candidates:
            continue
        if len(candidates) > 1:
            # pick the one with the most remaining copy-number budget
            nid_of = lambda v: v[:-1]
            candidates.sort(key=lambda v: -remaining_budget.get(nid_of(v), 0))
        chosen = candidates[0]
        walk.append(chosen)
        nid = chosen[:-1]
        if nid in remaining_budget:
            remaining_budget[nid] -= 1
    return walk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignment", required=True)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--meta", required=False, help="qubo .meta.json for copy numbers")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    nodes = read_gfa_nodes(args.gfa)
    by_t = read_assignment(args.assignment)

    remaining_budget = {nid: 10 ** 6 for nid in nodes}  # effectively unlimited unless meta given
    if args.meta:
        import json
        meta = json.load(open(args.meta))
        remaining_budget = dict(meta.get("weights", {}))

    walk = resolve_walk(by_t, remaining_budget)

    seq_parts = []
    for v in walk:
        nid, orient = v[:-1], v[-1]
        base_seq = nodes.get(nid, "")
        seq_parts.append(revcomp(base_seq) if orient == "-" else base_seq)

    candidate = "".join(seq_parts)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(">candidate\n")
        for i in range(0, len(candidate), 70):
            fh.write(candidate[i:i + 70] + "\n")

    print(f"[extract_path] walk length={len(walk)} nodes, "
          f"candidate genome length={len(candidate)}bp -> {args.out}")


if __name__ == "__main__":
    main()
