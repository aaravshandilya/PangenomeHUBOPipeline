#!/usr/bin/env python3
"""
Stage 2 - Read mapping: align short reads to the pangenome and tag nodes
with kmer counts / edge-traversal counts.

annotator = "kmer2node"    -> our own tool (real implementation below)
annotator = "minigraph"    -> real `minigraph -cx lr` binary if on PATH,
                              parsing its GAF path column; else falls back
                              to kmer2node
annotator = "graphaligner" -> real `GraphAligner` binary if on PATH,
                              parsing its GAF; else falls back to kmer2node

Output: an annotated GFA-like TSV: node_id, length, kmer_count, and an
edge-count TSV: from_id, to_id, traversal_count.

Usage:
    map_reads.py --gfa pangenome.gfa --fastq reads.fastq \
                  --annotator kmer2node --k 21 --out-prefix annotated
"""
import argparse
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


def read_gfa(path):
    nodes = {}
    edges = []
    for line in open(path):
        if line.startswith("S\t"):
            _, nid, seq = line.rstrip("\n").split("\t")[:3]
            nodes[nid] = seq
        elif line.startswith("L\t"):
            parts = line.rstrip("\n").split("\t")
            a, oa, b, ob = parts[1], parts[2], parts[3], parts[4]
            edges.append((a, oa, b, ob))
    return nodes, edges


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_fastq(path):
    reads = []
    with open(path) as fh:
        lines = fh.readlines()
    for i in range(0, len(lines), 4):
        if i + 1 >= len(lines):
            break
        reads.append(lines[i + 1].strip())
    return reads


def build_kmer_index(nodes, k):
    """Map each kmer -> list of (node_id, offset). Track uniqueness."""
    index = defaultdict(list)
    for nid, seq in nodes.items():
        for i in range(len(seq) - k + 1):
            index[seq[i:i + k]].append((nid, i))
    return index


def kmer2node_annotate(nodes, edges, reads, k=21):
    """
    Our tool: index node sequences by kmer, map each read's kmers (both
    strands) to nodes via unique-kmer voting, then tally per-node kmer
    hits and per-read best-node-order to approximate edge traversals.
    """
    index = build_kmer_index(nodes, k)
    node_kmer_hits = Counter()
    edge_hits = Counter()

    for seq in reads:
        for candidate in (seq, revcomp(seq)):
            if len(candidate) < k:
                continue
            hit_nodes_in_order = []
            for i in range(0, len(candidate) - k + 1, max(1, k // 2)):
                kmer = candidate[i:i + k]
                hits = index.get(kmer)
                if not hits:
                    continue
                # only count unique (non-repetitive) kmer hits towards a node,
                # matching the paper's discussion of unique-kmer expectation
                if len(hits) == 1:
                    nid, _ = hits[0]
                    node_kmer_hits[nid] += 1
                    if not hit_nodes_in_order or hit_nodes_in_order[-1] != nid:
                        hit_nodes_in_order.append(nid)
            for a, b in zip(hit_nodes_in_order, hit_nodes_in_order[1:]):
                edge_hits[(a, b)] += 1

    return node_kmer_hits, edge_hits


def try_real_annotator(binary_name, gfa_path, fastq_path, out_gaf, extra_args):
    if shutil.which(binary_name) is None:
        return False
    with open(out_gaf, "w") as out_fh:
        cmd = [binary_name, *extra_args, gfa_path, fastq_path]
        subprocess.run(cmd, stdout=out_fh, check=True)
    return True


def parse_gaf_to_counts(gaf_path, nodes):
    """GAF path column looks like '>12>5<3...' ; tally node visits + edges."""
    node_kmer_hits = Counter()
    edge_hits = Counter()
    node_re = re.compile(r"[><](\d+)")
    for line in open(gaf_path):
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 6:
            continue
        path_str = cols[5]
        visited = node_re.findall(path_str)
        for nid in visited:
            if nid in nodes:
                node_kmer_hits[nid] += 1
        for a, b in zip(visited, visited[1:]):
            edge_hits[(a, b)] += 1
    return node_kmer_hits, edge_hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--annotator", required=True,
                    choices=["kmer2node", "minigraph", "graphaligner"])
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    nodes, edges = read_gfa(args.gfa)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    node_hits, edge_hits = None, None

    if args.annotator == "minigraph":
        gaf = out_prefix.with_suffix(".gaf")
        if try_real_annotator("minigraph", args.gfa, args.fastq, gaf,
                               ["-cx", "lr"]):
            node_hits, edge_hits = parse_gaf_to_counts(gaf, nodes)
            print("[map_reads] annotated with real `minigraph`")
    elif args.annotator == "graphaligner":
        gaf = out_prefix.with_suffix(".gaf")
        if try_real_annotator("GraphAligner", args.gfa, args.fastq, gaf,
                               ["-g", args.gfa, "-f", args.fastq, "-a"]):
            node_hits, edge_hits = parse_gaf_to_counts(gaf, nodes)
            print("[map_reads] annotated with real `GraphAligner`")

    if node_hits is None:
        if args.annotator != "kmer2node":
            print(f"[map_reads] real `{args.annotator}` binary not found; "
                  f"falling back to kmer2node as a stand-in annotator")
        reads = read_fastq(args.fastq)
        node_hits, edge_hits = kmer2node_annotate(nodes, edges, reads, k=args.k)

    with open(f"{out_prefix}.nodes.tsv", "w") as fh:
        fh.write("node_id\tlength\tkmer_count\n")
        for nid, seq in nodes.items():
            fh.write(f"{nid}\t{len(seq)}\t{node_hits.get(nid, 0)}\n")

    with open(f"{out_prefix}.edges.tsv", "w") as fh:
        fh.write("from_id\tto_id\ttraversal_count\n")
        for (a, b), c in edge_hits.items():
            fh.write(f"{a}\t{b}\t{c}\n")

    print(f"[map_reads] annotator={args.annotator}: "
          f"{sum(1 for v in node_hits.values() if v > 0)}/{len(nodes)} nodes hit")


if __name__ == "__main__":
    main()
