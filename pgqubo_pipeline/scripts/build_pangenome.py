#!/usr/bin/env python3
"""
Stage 1b - Problem creation: build the pangenome graph (GFA).

Preferred path: shell out to the real `minigraph` binary (Li et al. 2020),
producing a broadly-colinear pangenome as used by the HPRC and by this paper
(Methods -> Pangenome creation).

Fallback path (no minigraph on PATH): a simplified progressive-alignment
graph builder. It threads each genome through the growing graph using
longest-common-substring anchors (difflib), reusing shared nodes and
creating new "bubble" nodes for novel/divergent sequence. This is NOT a
faithful reimplementation of minigraph -- it exists so the pipeline is
runnable end-to-end without the real binary -- but it produces a valid
GFA with the same semantics (S = node sequence, L = edge) that all
downstream stages consume.

Usage:
    build_pangenome.py --fasta pangenome_build_genomes.fasta --out pangenome.gfa
"""
import argparse
import difflib
import shutil
import subprocess
from pathlib import Path


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


def try_real_minigraph(fasta_path, out_gfa):
    if shutil.which("minigraph") is None:
        return False
    records = read_fasta(fasta_path)
    if len(records) < 1:
        return False
    per_genome_fastas = []
    tmp_dir = Path(out_gfa).parent / "_minigraph_tmp"
    tmp_dir.mkdir(exist_ok=True)
    for name, seq in records:
        p = tmp_dir / f"{name}.fa"
        with open(p, "w") as fh:
            fh.write(f">{name}\n{seq}\n")
        per_genome_fastas.append(str(p))
    cmd = ["minigraph", "-cxggs", *per_genome_fastas]
    with open(out_gfa, "w") as out_fh:
        subprocess.run(cmd, stdout=out_fh, check=True)
    return True


class GraphBuilder:
    """Minimal progressive pangenome builder (fallback, non-minigraph)."""

    def __init__(self, min_anchor=15):
        self.min_anchor = min_anchor
        self.nodes = {}      # node_id -> sequence
        self.edges = set()   # (from_id, from_orient, to_id, to_orient)
        self.next_id = 1
        self.backbone = []   # list of node ids representing genome 0's path

    def _new_node(self, seq):
        nid = self.next_id
        self.nodes[nid] = seq
        self.next_id += 1
        return nid

    def add_first_genome(self, seq, chunk=500):
        path = []
        for i in range(0, len(seq), chunk):
            nid = self._new_node(seq[i:i + chunk])
            path.append(nid)
        for a, b in zip(path, path[1:]):
            self.edges.add((a, "+", b, "+"))
        self.backbone = path
        return path

    def _backbone_seq_and_offsets(self):
        seq_parts, offsets = [], []
        pos = 0
        for nid in self.backbone:
            offsets.append((pos, nid))
            seq_parts.append(self.nodes[nid])
            pos += len(self.nodes[nid])
        return "".join(seq_parts), offsets

    def add_genome(self, seq):
        """Align `seq` against the current backbone using LCS anchors;
        reuse matching backbone nodes, insert new nodes for novel stretches."""
        backbone_seq, offsets = self._backbone_seq_and_offsets()
        sm = difflib.SequenceMatcher(None, backbone_seq, seq, autojunk=False)
        blocks = [b for b in sm.get_matching_blocks() if b.size >= self.min_anchor]

        path = []
        prev_b_end, prev_q_end = 0, 0
        for blk in blocks:
            # novel query sequence before this anchor -> new node(s)
            if blk.b > prev_q_end:
                novel = seq[prev_q_end:blk.b]
                if novel:
                    nid = self._new_node(novel)
                    path.append(nid)
            # reused backbone region -> find which backbone node(s) this spans
            span_start, span_end = blk.a, blk.a + blk.size
            for pos, nid in offsets:
                node_len = len(self.nodes[nid])
                if pos < span_end and pos + node_len > span_start:
                    path.append(nid)
            prev_b_end, prev_q_end = span_end, blk.b + blk.size

        if prev_q_end < len(seq):
            novel = seq[prev_q_end:]
            if novel:
                nid = self._new_node(novel)
                path.append(nid)

        # dedupe consecutive repeats of the same node id
        dedup = [path[0]] if path else []
        for nid in path[1:]:
            if nid != dedup[-1]:
                dedup.append(nid)
        path = dedup

        for a, b in zip(path, path[1:]):
            self.edges.add((a, "+", b, "+"))
        return path

    def write_gfa(self, out_path):
        with open(out_path, "w") as fh:
            fh.write("H\tVN:Z:1.0\n")
            for nid, seq in sorted(self.nodes.items()):
                fh.write(f"S\t{nid}\t{seq}\n")
            for a, oa, b, ob in sorted(self.edges):
                fh.write(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M\n")


def fallback_builder(fasta_path, out_gfa):
    records = read_fasta(fasta_path)
    if not records:
        raise ValueError(f"No sequences found in {fasta_path}")
    gb = GraphBuilder()
    gb.add_first_genome(records[0][1])
    for _, seq in records[1:]:
        gb.add_genome(seq)
    gb.write_gfa(out_gfa)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    used_real = try_real_minigraph(args.fasta, args.out)
    if used_real:
        print("[build_pangenome] built GFA with real `minigraph` binary")
    else:
        print("[build_pangenome] `minigraph` not found on PATH; "
              "using fallback progressive-alignment builder")
        fallback_builder(args.fasta, args.out)


if __name__ == "__main__":
    main()
