#!/usr/bin/env python3
"""
Stage 5b - Optional refinement: realign the short reads back to the
candidate sequence and take a consensus, correcting for point mutations
or small indels absorbed into the minigraph node-collapse step.

Preferred: real `minimap2`/`bwa` + `samtools mpileup` if on PATH.
Fallback: a simple k-mer-anchored pileup consensus (majority-vote base
at each anchored position; falls back to the candidate base if a
position isn't covered).

Usage:
    refine_consensus.py --candidate candidate.fasta --fastq reads.fastq \
                         --out candidate_opt.fasta
"""
import argparse
import shutil
import subprocess
from collections import Counter
from pathlib import Path


def read_fasta_single(path):
    seq = []
    for line in open(path):
        line = line.strip()
        if line and not line.startswith(">"):
            seq.append(line)
    return "".join(seq)


def read_fastq(path):
    reads = []
    lines = open(path).readlines()
    for i in range(0, len(lines), 4):
        if i + 1 >= len(lines):
            break
        reads.append(lines[i + 1].strip())
    return reads


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def try_real_pipeline(candidate_fa, fastq, out_fa):
    if not (shutil.which("minimap2") and shutil.which("samtools")):
        return False
    tmp = Path(out_fa).parent / "_realign_tmp"
    tmp.mkdir(exist_ok=True, parents=True)
    sam = tmp / "aln.sam"
    bam = tmp / "aln.sorted.bam"
    subprocess.run(["minimap2", "-a", candidate_fa, fastq], stdout=open(sam, "w"), check=True)
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    consensus_fa = tmp / "consensus.fa"
    subprocess.run(
        ["samtools", "consensus", "-f", "fasta", "-o", str(consensus_fa), str(bam)],
        check=True,
    )
    shutil.copy(consensus_fa, out_fa)
    return True


def kmer_anchored_consensus(candidate, reads, k=21):
    """Lightweight fallback consensus: anchor each read to the candidate via
    unique k-mer matches, then majority-vote the base at each covered
    position."""
    index = {}
    for i in range(len(candidate) - k + 1):
        kmer = candidate[i:i + k]
        index.setdefault(kmer, []).append(i)

    votes = [Counter() for _ in range(len(candidate))]
    for read in reads:
        for strand in (read, revcomp(read)):
            if len(strand) < k:
                continue
            anchor_offsets = []
            for i in range(0, len(strand) - k + 1, max(1, k // 2)):
                kmer = strand[i:i + k]
                hits = index.get(kmer)
                if hits and len(hits) == 1:
                    anchor_offsets.append((i, hits[0]))
            if not anchor_offsets:
                continue
            # use the first anchor to place the whole read (ungapped)
            read_off, cand_off = anchor_offsets[0]
            start_in_cand = cand_off - read_off
            for j, base in enumerate(strand):
                pos = start_in_cand + j
                if 0 <= pos < len(candidate):
                    votes[pos][base] += 1

    consensus = []
    for pos, base in enumerate(candidate):
        if votes[pos]:
            consensus.append(votes[pos].most_common(1)[0][0])
        else:
            consensus.append(base)
    return "".join(consensus)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    if try_real_pipeline(args.candidate, args.fastq, args.out):
        print("[refine_consensus] used real minimap2 + samtools consensus")
        return

    print("[refine_consensus] minimap2/samtools not found; "
          "using fallback k-mer-anchored consensus")
    candidate = read_fasta_single(args.candidate)
    reads = read_fastq(args.fastq)
    consensus = kmer_anchored_consensus(candidate, reads)

    with open(args.out, "w") as fh:
        fh.write(">candidate_opt\n")
        for i in range(0, len(consensus), 70):
            fh.write(consensus[i:i + 70] + "\n")


if __name__ == "__main__":
    main()
