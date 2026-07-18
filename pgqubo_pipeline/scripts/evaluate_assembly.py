#!/usr/bin/env python3
"""
Stage 5c - Classical postprocessing: evaluate a candidate assembly against
ground truth, reproducing the paper's Table 1/2/3 columns:
    %Covered, %Used, Contigs, Breaks, Indels, No.Diff, N50, %Identity

Preferred: real `bwa mem` (candidate vs truth), parsing primary +
supplementary alignments as "contigs", exactly as described in Methods
("Classical postprocessing: evaluating path solutions").

Fallback: difflib.SequenceMatcher matching blocks between candidate and
truth stand in for alignment blocks/contigs; gaps >= break_min_gap in
either sequence split a contig (a "break"); mismatches within a block are
substitutions; short unmatched runs within a block are indels or
diff-regions depending on length, using the config thresholds.

Usage:
    evaluate_assembly.py --candidate candidate.fasta --truth test_genome.fasta \
                          --config config.yaml --out metrics.json
"""
import argparse
import difflib
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import yaml


def read_fasta_single(path):
    seq = []
    for line in open(path):
        line = line.strip()
        if line and not line.startswith(">"):
            seq.append(line)
    return "".join(seq)


def n50(lengths):
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    running = 0
    for length in lengths:
        running += length
        if running >= total / 2:
            return length
    return lengths[-1]


def try_real_bwa(candidate_fa, truth_fa, cfg):
    if not (shutil.which("bwa") and shutil.which("samtools")):
        return None
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        shutil.copy(truth_fa, tmp / "truth.fa")
        subprocess.run(["bwa", "index", str(tmp / "truth.fa")], check=True,
                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sam = tmp / "aln.sam"
        subprocess.run(["bwa", "mem", str(tmp / "truth.fa"), candidate_fa],
                        stdout=open(sam, "w"), stderr=subprocess.DEVNULL, check=True)
        # NOTE: parsing CIGAR/primary+supplementary alignments into the
        # exact metric set is left as an integration point; a full SAM
        # parser (pysam) is the natural way to complete this branch.
        return None  # fall through to the difflib approximation for now
    return None


def evaluate_with_difflib(candidate, truth, cfg):
    ecfg = cfg["evaluation"]
    sm = difflib.SequenceMatcher(None, truth, candidate, autojunk=False)
    blocks = [b for b in sm.get_matching_blocks() if b.size > 0]

    # merge blocks into "contigs": consecutive blocks are part of the same
    # contig unless the gap in either sequence exceeds break_min_gap
    contigs = []
    current = [blocks[0]] if blocks else []
    breaks = 0
    for prev, cur in zip(blocks, blocks[1:]):
        gap_truth = cur.a - (prev.a + prev.size)
        gap_cand = cur.b - (prev.b + prev.size)
        if gap_truth > ecfg["break_min_gap"] or gap_cand > ecfg["break_min_gap"]:
            contigs.append(current)
            current = [cur]
            breaks += 1
        else:
            current.append(cur)
    if current:
        contigs.append(current)

    contig_lengths = []
    covered_bp = 0
    used_bp = 0
    indels = 0
    diffs = 0
    matched_bp = 0
    total_block_span = 0

    for group in contigs:
        span_len = (group[-1].a + group[-1].size) - group[0].a
        contig_lengths.append(max(span_len, 1))
        for blk in group:
            covered_bp += blk.size
            used_bp += blk.size
            matched_bp += blk.size
        for prev, cur in zip(group, group[1:]):
            gap = cur.b - (prev.b + prev.size)
            if gap >= ecfg["indel_min_len"]:
                indels += 1
            elif gap > 0:
                diffs += 1

    pct_covered = 100.0 * covered_bp / max(len(truth), 1)
    pct_used = 100.0 * used_bp / max(len(candidate), 1)
    total_block_span = sum(contig_lengths)
    pct_identity = 100.0 * matched_bp / max(total_block_span, 1)

    return {
        "pct_covered": round(pct_covered, 2),
        "pct_used": round(pct_used, 2),
        "n_contigs": len(contigs),
        "breaks": breaks,
        "indels": indels,
        "n_diff": diffs,
        "n50": n50(contig_lengths),
        "pct_identity": round(pct_identity, 2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="", help="e.g. seed=1;annotator=kmer2node;solver=gurobi")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))
    candidate = read_fasta_single(args.candidate)
    truth = read_fasta_single(args.truth)

    metrics = try_real_bwa(args.candidate, args.truth, cfg)
    if metrics is None:
        metrics = evaluate_with_difflib(candidate, truth, cfg)
        metrics["method"] = "difflib_fallback"
    else:
        metrics["method"] = "bwa_mem"

    metrics["tag"] = args.tag
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    json.dump(metrics, open(args.out, "w"), indent=2)
    print(f"[evaluate_assembly] {args.tag}: {metrics}")


if __name__ == "__main__":
    main()
