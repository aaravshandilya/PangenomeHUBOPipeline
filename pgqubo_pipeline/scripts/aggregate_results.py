#!/usr/bin/env python3
"""
Stage 5e - Aggregate all per-instance metrics.json files into:
    results/summary_by_annotator_solver.csv  (mirrors Table 2/3: mean (std))
    results/summary.png                       (bar comparison across solvers)

Usage:
    aggregate_results.py --metrics-glob "results/**/metrics.json" --out-dir results
"""
import argparse
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_tag(tag):
    out = {}
    for kv in tag.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics-glob", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    rows = []
    for path in glob.glob(args.metrics_glob, recursive=True):
        m = json.load(open(path))
        row = parse_tag(m.get("tag", ""))
        row.update({k: v for k, v in m.items() if k not in ("tag",)})
        rows.append(row)

    if not rows:
        print("[aggregate_results] no metrics files found")
        return

    df = pd.DataFrame(rows)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "all_instances.csv", index=False)

    numeric_cols = ["pct_covered", "pct_used", "n_contigs", "breaks",
                    "indels", "n_diff", "n50", "pct_identity"]
    numeric_cols = [c for c in numeric_cols if c in df.columns]

    group_cols = [c for c in ("annotator", "solver") if c in df.columns]
    if group_cols:
        summary = df.groupby(group_cols)[numeric_cols].agg(["mean", "std"])
        summary.to_csv(out_dir / "summary_by_annotator_solver.csv")
        print(summary)

        # simple bar chart: mean n_contigs and mean pct_identity per group
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        means = df.groupby(group_cols)[numeric_cols].mean()
        means["n_contigs"].plot(kind="bar", ax=axes[0], title="Mean contigs")
        means["pct_identity"].plot(kind="bar", ax=axes[1], title="Mean %identity")
        plt.tight_layout()
        plt.savefig(out_dir / "summary.png", dpi=150)
        print(f"[aggregate_results] wrote {out_dir/'summary.png'} and "
              f"{out_dir/'summary_by_annotator_solver.csv'}")


if __name__ == "__main__":
    main()
