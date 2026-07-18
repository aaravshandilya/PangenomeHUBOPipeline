# Pangenome-Guided Sequence Assembly — Snakemake Pipeline

A Snakemake reimplementation of the 5-stage workflow from
**Strelchuk et al., "Pangenome-guided sequence assembly via binary optimization,"
Briefings in Bioinformatics 2026 (bbag084)**.

```
1. Problem creation       generate_population.py, build_pangenome.py, shred_genome.py
2. Read mapping           map_reads.py           (kmer2node / minigraph / GraphAligner)
3. Copy number estimation estimate_copy_number.py (pathfinder or fallback)
4. Path finding            build_qubo.py + solve_qubo.py
                           (classical_sa / gurobi / mqlib / pathfinder / dwave)
5. Solution processing     extract_path.py, refine_consensus.py, evaluate_assembly.py
                           -> aggregate_results.py (Table 1/2/3-style summaries)
```

## Quick start

```bash
conda env create -f envs/environment.yml
conda activate pgqubo
snakemake -j 8 --configfile config.yaml -n     # dry run, inspect the DAG
snakemake -j 8 --configfile config.yaml         # run everything
```

Edit `config.yaml` to set the seed/replicate count, which annotators and
solvers to sweep, and the population/sequencing/QUBO parameters (Lagrange
multipliers λ1, λ2, walk-budget α, etc., matching the paper's Sec.
"Mapping guided alignment to binary optimization").

Outputs land in `results/`:
- `results/seed_<n>/pangenome.gfa` — the synthetic pangenome graph
- `results/seed_<n>/annot_<tool>/test_<i>/...` — per read-mapper artifacts
- `.../solver_<solver>_tl<t>/candidate_opt.fasta` — final assembled sequence
- `.../metrics.json` — %Covered, %Used, Contigs, Breaks, Indels, N50, %Identity
- `results/summary_by_annotator_solver.csv` and `summary.png` — aggregated,
  paper-Table-2-style comparison across every annotator × solver combination

## What's real vs. what's a stand-in

This pipeline is designed to be **runnable out of the box** even without the
paper's specialized (and non-pip-installable) dependencies, while giving you
a clean drop-in point to swap in the real tools:

| Stage | Real tool (used automatically if on `$PATH` / installed) | Built-in fallback (always works) |
|---|---|---|
| Pangenome construction | `minigraph` | simplified progressive-alignment graph builder (`build_pangenome.py`) |
| Read mapping | `minigraph`, `GraphAligner` | our own `kmer2node` k-mer index/voting implementation (real, not a stub) |
| Copy number estimation | Oatk's `pathfinder` | modal-depth normalization |
| QUBO solving | `gurobipy` (needs license), D-Wave `dwave-ocean-sdk` (needs Leap token), MQLib binary | classical simulated annealing over the exact same sparse QUBO (`solve_qubo.py`) |
| Realignment/polish | `minimap2` + `samtools consensus` | k-mer-anchored majority-vote consensus |
| Evaluation | `bwa mem` (hook present, CIGAR parsing left as an integration point — see comment in `evaluate_assembly.py`) | `difflib`-based alignment-block approximation of contigs/breaks/indels/identity |

The **QUBO formulation itself is a full, real implementation** of the
paper's equations (2), (4)–(6): oriented tangle resolution with a virtual
"end" node, `λ1` one-vertex-per-timestep constraint, `λ2` valid-edge-traversal
constraint, and the copy-number-matching cost term, all built as a sparse
`{(i,j): coeff}` dictionary so any of Gurobi/MQLib/D-Wave/simulated annealing
can consume the identical problem instance (`build_qubo.py`, `solve_qubo.py`).

To wire in the real bioinformatics binaries, just install them
(`conda install -c bioconda minigraph graphaligner`, build Oatk's
`pathfinder`, `pip install gurobipy` with a license, `pip install
dwave-ocean-sdk` with a Leap API token) — each script auto-detects the
binary/library and only falls back if it's missing, so no config changes
are required.

## Notes on scale

Table 2 of the paper uses pangenomes with ~45–80 nodes and QUBOs with
`(2N+1)·T` binary variables. At that scale, exact classical solvers
(Gurobi/MQLib) and hybrid quantum annealers are appropriate; the bundled
`classical_sa` fallback will still run but is slower and lower-quality
than a tuned tabu search — treat it as a correctness baseline, not a
performance comparison, unless you install MQLib/Gurobi.

## Extending

- **Diploid tangle resolution** (paper's Problem 3): duplicate the variable
  set in `build_qubo.py` as `y_{t,v}` sharing the `C3` term — noted as a
  `TODO`-style extension point in that script's docstring.
- **Edge-weighted QUBO**: the paper notes (Discussion) that ignoring edge
  weights costs extra "breaks" from unresolved inversions; `build_qubo.py`'s
  `C2` term is the natural place to add an edge-weight-dependent bonus.
