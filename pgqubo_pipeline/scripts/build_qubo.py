#!/usr/bin/env python3
"""
Stage 4a - Formulate oriented tangle resolution as a QUBO (paper eqs 2, 4-6).

Variables: x_{t,v} for t = 1..T, v in V+ ∪ V- ∪ {end}.
  T = ceil(alpha * sum_v w(v))                     (walk-length budget)

Cost = lambda1 * C1 (exactly one node per time step)
     + lambda2 * C2 (only traverse real graph edges, or stay in `end`)
     + C3           (visit counts match copy numbers, eq. 2)

Output: a sparse QUBO as (i, j, coeff) triplets plus a JSON variable map,
so any solver (classical SA, Gurobi, MQLib, D-Wave/dimod, QAOA) can consume
the exact same problem instance.

Usage:
    build_qubo.py --gfa pangenome.gfa --copy-numbers copy_numbers.tsv \
                   --config config.yaml --out-prefix qubo/instance
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import yaml


def read_gfa_edges(path):
    edges = set()
    node_ids = []
    for line in open(path):
        if line.startswith("S\t"):
            node_ids.append(line.split("\t")[1])
        elif line.startswith("L\t"):
            parts = line.rstrip("\n").split("\t")
            a, oa, b, ob = parts[1], parts[2], parts[3], parts[4]
            edges.add((a, oa, b, ob))
    return node_ids, edges


def read_copy_numbers(path):
    w = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            nid, _depth, cn = line.strip().split("\t")
            w[nid] = int(cn)
    return w


def oriented_vertices(node_ids):
    """v+ and v- for each unsigned node id."""
    return [f"{nid}+" for nid in node_ids] + [f"{nid}-" for nid in node_ids]


def oriented_edge_set(edges):
    """Expand undirected/oriented GFA edges (a,oa,b,ob) into directed pairs
    over the signed vertex set, plus their reverse-complement pair, matching
    the paper's convention that an edge A+->B+ implies B-->A-."""
    signed_edges = set()
    for a, oa, b, ob in edges:
        sa = f"{a}{oa}"
        sb = f"{b}{ob}"
        signed_edges.add((sa, sb))
        # reverse complement implied edge
        flip = lambda s: s[:-1] + ("-" if s[-1] == "+" else "+")
        signed_edges.add((flip(sb), flip(sa)))
    return signed_edges


class VarIndex:
    """Maps (t, v) -> integer variable index. v ranges over signed vertices
    plus the special string 'end'."""

    def __init__(self, T, vertices):
        self.T = T
        self.vertices = vertices + ["end"]
        self._idx = {}
        i = 0
        for t in range(1, T + 1):
            for v in self.vertices:
                self._idx[(t, v)] = i
                i += 1
        self.n_vars = i

    def __getitem__(self, key):
        return self._idx[key]


def add_q(Q, i, j, val):
    if i > j:
        i, j = j, i
    Q[(i, j)] += val


def build_qubo(node_ids, edges, weights, alpha, lambda1, lambda2):
    vertices = oriented_vertices(node_ids)
    signed_edges = oriented_edge_set(edges)

    total_w = sum(weights.get(nid, 0) for nid in node_ids)
    T = max(1, int(alpha * max(total_w, 1)))

    idx = VarIndex(T, vertices)
    Q = defaultdict(float)

    # --- C1: exactly one vertex (incl. `end`) visited at each time step ----
    # lambda1 * sum_t ( sum_v x_t,v + x_t,end - 1 )^2
    for t in range(1, T + 1):
        all_v = vertices + ["end"]
        var_ids = [idx[(t, v)] for v in all_v]
        # expand (sum x_i - 1)^2 = sum x_i (linear via x_i^2=x_i)
        #   + 2*sum_{i<j} x_i x_j  - 2*sum x_i + 1 (constant dropped)
        for vi in var_ids:
            add_q(Q, vi, vi, lambda1 * (1 - 2))  # diag: x_i^2 term + (-2 x_i) term folded (x_i^2=x_i)
        for a_pos in range(len(var_ids)):
            for b_pos in range(a_pos + 1, len(var_ids)):
                add_q(Q, var_ids[a_pos], var_ids[b_pos], lambda1 * 2)

    # --- C2: only traverse real edges between consecutive time steps -------
    # lambda2 * sum_{t=1}^{T-1} [ sum_{v,v' not edge} x_t,v x_{t+1,v'}
    #                              + sum_v x_t,end x_{t+1,v} ]
    for t in range(1, T):
        for va in vertices:
            for vb in vertices:
                if (va, vb) not in signed_edges:
                    i, j = idx[(t, va)], idx[(t + 1, vb)]
                    add_q(Q, i, j, lambda2)
        for vb in vertices:
            i, j = idx[(t, "end")], idx[(t + 1, vb)]
            add_q(Q, i, j, lambda2)

    # --- C3: visit counts match copy numbers (oriented: v+ and v- share w(v)) -
    # sum_v ( sum_t [x_t,v+ + x_t,v-] - w(v) )^2
    for nid in node_ids:
        wv = weights.get(nid, 0)
        vplus, vminus = f"{nid}+", f"{nid}-"
        var_ids = [idx[(t, vplus)] for t in range(1, T + 1)] + \
                  [idx[(t, vminus)] for t in range(1, T + 1)]
        for vi in var_ids:
            add_q(Q, vi, vi, (1 - 2 * wv))
        for a_pos in range(len(var_ids)):
            for b_pos in range(a_pos + 1, len(var_ids)):
                add_q(Q, var_ids[a_pos], var_ids[b_pos], 2)
        # constant term wv^2 dropped (doesn't affect argmin)

    return Q, idx, T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--copy-numbers", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))["qubo"]
    node_ids, edges = read_gfa_edges(args.gfa)
    weights = read_copy_numbers(args.copy_numbers)

    Q, idx, T = build_qubo(node_ids, edges, weights,
                            alpha=cfg["alpha"],
                            lambda1=cfg["lambda1"],
                            lambda2=cfg["lambda2"])

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    with open(f"{out_prefix}.qubo.tsv", "w") as fh:
        fh.write("i\tj\tcoeff\n")
        for (i, j), c in Q.items():
            if c != 0:
                fh.write(f"{i}\t{j}\t{c}\n")

    meta = {
        "n_vars": idx.n_vars,
        "T": T,
        "vertices": idx.vertices,
        "node_ids": node_ids,
        "weights": weights,
    }
    json.dump(meta, open(f"{out_prefix}.meta.json", "w"), indent=2)
    print(f"[build_qubo] {idx.n_vars} binary variables (T={T}, "
          f"{len(node_ids)} nodes x2 orientations + end), "
          f"{sum(1 for c in Q.values() if c != 0)} nonzero QUBO terms")


if __name__ == "__main__":
    main()
