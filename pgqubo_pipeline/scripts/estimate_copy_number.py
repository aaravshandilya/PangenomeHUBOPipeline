#!/usr/bin/env python3
"""
Stage 3 - Copy number estimation.

Preferred: shell out to the real `pathfinder` binary (from Oatk), which the
paper uses for copy-number estimation regardless of which annotator produced
the counts (see Discussion: "For copy number estimation, we only use the
capabilities of pathfinder").

Fallback: depth-normalize kmer counts by the modal per-base coverage
(a robust single-copy depth estimate), matching the "depth-normalized
counts" pathfinder would otherwise supply.

Usage:
    estimate_copy_number.py --nodes annotated.nodes.tsv --config config.yaml \
                             --out copy_numbers.tsv
"""
import argparse
import shutil
import subprocess
from pathlib import Path

import numpy as np
import yaml


def try_real_pathfinder(binary, nodes_tsv, out_tsv):
    if shutil.which(binary) is None:
        return False
    # Real Oatk pathfinder takes a GFA with depth tags; here we assume a
    # wrapper convention: `pathfinder estimate-cn <nodes_tsv> -o <out_tsv>`.
    # Adjust to match the actual CLI of your pathfinder build.
    cmd = [binary, "estimate-cn", nodes_tsv, "-o", out_tsv]
    subprocess.run(cmd, check=True)
    return True


def fallback_estimate(nodes_tsv, out_tsv):
    node_ids, lengths, counts = [], [], []
    with open(nodes_tsv) as fh:
        next(fh)  # header
        for line in fh:
            nid, length, kc = line.strip().split("\t")
            node_ids.append(nid)
            lengths.append(int(length))
            counts.append(int(kc))

    lengths = np.array(lengths, dtype=float)
    counts = np.array(counts, dtype=float)
    depth = np.divide(counts, np.maximum(lengths, 1))

    nonzero_depth = depth[depth > 0]
    if len(nonzero_depth) == 0:
        single_copy_depth = 1.0
    else:
        # robust single-copy depth estimate: median of the lower half of
        # nonzero depths (majority of nodes in a well-formed pangenome are
        # single copy relative to any one sample)
        single_copy_depth = max(np.median(nonzero_depth), 1e-6)

    copy_numbers = np.round(depth / single_copy_depth).astype(int)
    copy_numbers = np.clip(copy_numbers, 0, None)

    with open(out_tsv, "w") as fh:
        fh.write("node_id\tdepth\tcopy_number\n")
        for nid, d, cn in zip(node_ids, depth, copy_numbers):
            fh.write(f"{nid}\t{d:.4f}\t{cn}\n")

    print(f"[estimate_copy_number] fallback depth-normalization: "
          f"single_copy_depth~={single_copy_depth:.3f}, "
          f"{int((copy_numbers > 0).sum())}/{len(copy_numbers)} nodes with cn>0")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    binary = cfg["solver_params"]["pathfinder"]["binary"]
    if try_real_pathfinder(binary, args.nodes, args.out):
        print(f"[estimate_copy_number] used real `{binary}`")
    else:
        print(f"[estimate_copy_number] `{binary}` not found on PATH; using fallback")
        fallback_estimate(args.nodes, args.out)


if __name__ == "__main__":
    main()
