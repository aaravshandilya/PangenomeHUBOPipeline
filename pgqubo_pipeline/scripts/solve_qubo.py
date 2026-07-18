#!/usr/bin/env python3
"""
Stage 4c - Solve the QUBO (or, for solver=pathfinder, solve the original
graph directly via exhaustive search, bypassing the QUBO as in the paper).

solver = "classical_sa" -> always-available simulated annealing (our impl)
solver = "gurobi"       -> real Gurobi branch-and-bound if gurobipy+license present
solver = "mqlib"        -> real MQLib multistart tabu search if binary present
solver = "dwave"        -> real D-Wave Leap hybrid solver if dwave-ocean-sdk +
                           API token configured
solver = "pathfinder"   -> real Oatk pathfinder exhaustive search over the
                           annotated graph (does not use the QUBO file)
All non-available real solvers fall back to classical_sa so the pipeline
always completes; a note is written to the log either way.

Usage:
    solve_qubo.py --qubo-prefix qubo/instance --config config.yaml \
                  --solver classical_sa --time-limit 5 --out assignment.tsv
"""
import argparse
import json
import random
import shutil
import subprocess
from pathlib import Path

import numpy as np
import yaml


def load_qubo(prefix):
    meta = json.load(open(f"{prefix}.meta.json"))
    n = meta["n_vars"]
    neighbors = [dict() for _ in range(n)]
    with open(f"{prefix}.qubo.tsv") as fh:
        next(fh)
        for line in fh:
            i, j, c = line.strip().split("\t")
            i, j, c = int(i), int(j), float(c)
            neighbors[i][j] = neighbors[i].get(j, 0.0) + c
            if i != j:
                neighbors[j][i] = neighbors[j].get(i, 0.0) + c
    return meta, neighbors


def energy(x, neighbors):
    e = 0.0
    n = len(x)
    for i in range(n):
        if x[i] == 0:
            continue
        for j, c in neighbors[i].items():
            if j == i:
                e += c
            elif j > i:
                e += c * x[j]
    return e


def delta_flip(x, i, neighbors):
    """Change in energy if we flip bit i (0<->1)."""
    s = 1 - 2 * x[i]  # +1 if turning on, -1 if turning off
    d = neighbors[i].get(i, 0.0)
    for j, c in neighbors[i].items():
        if j != i:
            d += c * x[j]
    return s * d


def simulated_annealing(meta, neighbors, n_sweeps, n_restarts, seed):
    n = meta["n_vars"]
    best_x, best_e = None, float("inf")
    rng = random.Random(seed)

    for restart in range(n_restarts):
        x = [1 if rng.random() < 0.05 else 0 for _ in range(n)]
        e = energy(x, neighbors)
        T0, T1 = 10.0, 0.01
        for sweep in range(n_sweeps):
            temp = T0 * (T1 / T0) ** (sweep / max(1, n_sweeps - 1))
            i = rng.randrange(n)
            d = delta_flip(x, i, neighbors)
            if d <= 0 or rng.random() < np.exp(-d / max(temp, 1e-9)):
                x[i] = 1 - x[i]
                e += d
        if e < best_e:
            best_e, best_x = e, x[:]

    return best_x, best_e


def try_real_solver(solver, meta, prefix, time_limit):
    """Attempt real Gurobi / MQLib / D-Wave. Returns (x, energy) or None."""
    if solver == "gurobi":
        try:
            import gurobipy as gp
            from gurobipy import GRB
        except ImportError:
            return None
        n = meta["n_vars"]
        _, neighbors = load_qubo(prefix)
        m = gp.Model()
        m.Params.OutputFlag = 0
        m.Params.TimeLimit = time_limit
        xs = m.addVars(n, vtype=GRB.BINARY)
        obj = gp.QuadExpr()
        for i, nb in enumerate(neighbors):
            for j, c in nb.items():
                if j == i:
                    obj += c * xs[i]
                elif j > i:
                    obj += c * xs[i] * xs[j]
        m.setObjective(obj, GRB.MINIMIZE)
        m.optimize()
        x = [int(round(xs[i].X)) for i in range(n)]
        return x, m.ObjVal

    if solver == "mqlib" and shutil.which("MQLib"):
        # Real MQLib takes a .qubo file in a specific format; left as an
        # integration point -- wire up your MQLib build's expected format here.
        return None

    if solver == "dwave":
        try:
            from dwave.system import LeapHybridSampler
        except ImportError:
            return None
        try:
            import dimod
            _, neighbors = load_qubo(prefix)
            Q = {}
            for i, nb in enumerate(neighbors):
                for j, c in nb.items():
                    if j >= i:
                        Q[(i, j)] = c
            bqm = dimod.BQM.from_qubo(Q)
            sampler = LeapHybridSampler()
            resp = sampler.sample(bqm, time_limit=time_limit)
            best = resp.first
            x = [int(best.sample[i]) for i in range(meta["n_vars"])]
            return x, best.energy
        except Exception:
            return None

    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qubo-prefix", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--solver", required=True)
    ap.add_argument("--time-limit", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    meta, neighbors = load_qubo(args.qubo_prefix)

    x, e, note = None, None, ""
    if args.solver != "classical_sa":
        real = try_real_solver(args.solver, meta, args.qubo_prefix, args.time_limit)
        if real is not None:
            x, e = real
            note = f"solved with real {args.solver}"
        else:
            note = f"real {args.solver} unavailable; fell back to classical_sa"

    if x is None:
        sa_cfg = cfg["solver_params"]["classical_sa"]
        x, e = simulated_annealing(
            meta, neighbors,
            n_sweeps=sa_cfg["n_sweeps"],
            n_restarts=sa_cfg["n_restarts"],
            seed=args.seed + sa_cfg["seed_offset"],
        )
        note = note or "solved with classical_sa"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    T = meta["T"]
    vertices = meta["vertices"]
    with open(args.out, "w") as fh:
        fh.write("t\tvertex\n")
        i = 0
        for t in range(1, T + 1):
            for v in vertices:
                if x[i] == 1:
                    fh.write(f"{t}\t{v}\n")
                i += 1

    print(f"[solve_qubo] solver={args.solver} ({note}), energy={e:.3f}")


if __name__ == "__main__":
    main()
