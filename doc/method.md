# Method

This document describes the algorithmic method used by Sharda to produce phased
haplotype sequences from aligned short reads. The same backbone/DBG/unitig
pipeline also supports an SV-oriented execution mode that calls simple indels
from the unitig graph instead of running haplotype flow decomposition.

The tandem repeat BED is optional. When no TR annotations are supplied, the
pipeline still runs, but all reads are effectively treated as ORR reads and the
IRR anchor-chaining path is never activated.

## Overview

Sharda currently has three execution modes over the same graph-construction
pipeline:

- **Default combined mode** — build the graph, compact to unitigs, call
   simple indels from the unitig graph, run ILP flow decomposition, and emit
   both haplotype FASTA and SV outputs.
- **SV-only mode** (`--sv-only`) — build the graph, compact to unitigs, call
   simple indels from the unitig graph, emit the SV outputs, and skip
   haplotype flow decomposition.
- **Haplotype-only mode** (`--hap-only`) — build the graph, compact to
   unitigs, run ILP flow decomposition, and skip SV calling plus SV-oriented
   unitig outputs.
- **Unitig-only mode** (`--unitig-only`) — stop after unitig graph
  construction and artifact emission. This mode skips both haplotype flow
  decomposition and SV calling.

Sharda follows this pipeline for each target region:

1. **Backbone construction** — build a linear chain of positional de Bruijn
   graph nodes from the reference sequence.
2. **Read addition** — thread aligned reads through the graph, choosing between
   positional (ORR) and anchor-chained (IRR) strategies.
3. **Graph cleaning** — iteratively remove under-supported tips and low-weight edges; bubble popping is currently disabled.
4. **Unitig compaction** — collapse maximal non-branching paths into unitigs.
5. **SV calling and/or flow decomposition** — depending on the execution mode,
   call simple indels from alternate unitig traversals and/or solve an ILP to
   decompose unitig coverage into haplotype paths.
6. **Output** — emit the outputs selected by the mode: haplotype FASTA,
   `<prefix>.sv.vcf`, SV-oriented unitig artifacts, graph artifacts, or only
   unitig/debug artifacts.

---

## 1. Positional de Bruijn graph

A standard de Bruijn graph maps each k-mer to a single node, which is
problematic in tandem repeat regions where the same k-mer appears at multiple
genomic positions. Sharda's graph is *positional*: backbone nodes are keyed by
their reference coordinate, while non-backbone (read-derived) nodes are keyed
by k-mer string and can accumulate multiple implied coordinates in a coordinate
set. This preserves the linear order of the reference while still allowing one
read-derived node to represent the same novel k-mer reused at multiple nearby
placements.

Each node stores:

| Field | Meaning |
| ----- | ------- |
| `id` | Unique 64-bit identifier |
| `kmer` | The k-mer string (length *k*) |
| `ref_pos` | Primary local coordinate; exact for backbone, minimum implied coordinate for non-backbone nodes |
| `ref_positions` | All stored local coordinates associated with the node |
| `is_backbone` | Whether the node was derived from the reference |
| `depth` | Number of reads covering this node |

Edges are directed, weighted by the number of reads supporting the transition.

## 2. Backbone construction

The reference sequence for the target region is decomposed into overlapping
k-mers at consecutive positions. Each k-mer becomes a backbone node with a
unique `ref_pos`, and edges connect consecutive positions:

```text
ref[0..k-1]  →  ref[1..k]  →  ref[2..k+1]  →  ...
```

Tandem repeat annotations (from the TR BED file) are mapped to the set of
backbone nodes whose positions fall within each repeat, enabling per-TR anchor
lookup during read addition.

## 3. Read classification

Before adding a read to the graph, it is classified based on alignment evidence:

### Evidence detection

A read is an **evidence read** if it has any of:

- **Soft-clipping** at either end (≥ 5 bp) indicating the read extends beyond
  the reference alignment.
- **Large indels** (≥ 5 bp) in the CIGAR string.
- **Split alignment** (SA supplementary alignment tag).
- **Improper pair** flag set.
- **Tandem repeat overlap** — the aligned region overlaps an annotated TR.

### ORR vs IRR

| Type | Full name | Description |
| ---- | --------- | ----------- |
| **ORR** | Ordinary Reference Read | Aligns entirely within the reference backbone. Added using positional CIGAR walk. |
| **IRR** | In-Repeat Read | Overlaps a tandem repeat region. Added using anchor chaining. |

A read is classified as IRR when it is an evidence read *and* its aligned
extent overlaps a tandem repeat. All other reads (including non-evidence reads)
are treated as ORR. If no TR BED is provided, there are no TR overlaps, so all
reads remain ORR.

## 4. Read addition

### ORR path (implied-coordinate placement)

For an ORR read starting at local coordinate $N$, the first read k-mer is given
implied coordinate $N$ and each subsequent read k-mer gets implied coordinate
$N + i$.

At each read k-mer position:

1. Search backbone nodes carrying the same k-mer.
2. If one or more backbone nodes match, reuse the backbone occurrence whose
   coordinate is closest to the implied coordinate.
3. If no backbone node matches, reuse or create a non-backbone read node keyed
   by k-mer string and record the implied coordinate in that node's
   `ref_positions` set.
4. Add an edge from the previous chosen node. If the first read k-mer is novel
   and $N > 0$, also add the requested branch edge from backbone coordinate
   $N - 1$ into that first novel node.

This produces a path that follows the alignment start coordinate, but it is not
restricted to exact positional backbone matches when a repeated backbone k-mer
has a closer occurrence elsewhere.

### IRR path (anchor chaining)

For reads that overlap a tandem repeat, exact ORR-style implied placement can
still be ambiguous. Instead:

1. **Extract k-mers** from the read sequence.
2. **Find anchors** — k-mers that appear *exactly once* among the backbone
   nodes belonging to the overlapping TR. Uniqueness within the TR ensures an
   unambiguous mapping.
3. **Chain anchors** via O(n²) dynamic programming. The DP maximises the number
   of co-linear anchors (read position and backbone position both increasing)
   while penalising gap discrepancy:

   ```text
   score(i) = max_{j < i, bb_pos_j < bb_pos_i}  score(j) + 1 − |read_gap − bb_gap|
   ```

4. **Walk the read**: at each k-mer position, if an anchor is available, use
   the backbone node; otherwise, create a hash-based read node. This interleaves
   positional and standard DBG strategies.
5. **Fallback**: if no anchors are found, fall back to the ORR implied-coordinate path.

### Haplotype edges

After placing both reads of a pair, if either read is an evidence read, a
**haplotype edge** is added between the last node of the left read and the first
node of the right read (ordered by alignment position). Haplotype edges are not
ordinary graph edges — they do not participate in path traversal during graph
cleaning. Instead, they are used as constraints during flow decomposition to
enforce phasing.

## 5. Graph cleaning

Graph cleaning removes noise and errors iteratively (up to 10 rounds, stopping
when no changes are made). The current implementation performs tip removal and
low-weight edge pruning each round, and explicitly skips bubble popping.

### Tip removal

A **tip** is a dead-end path (in-degree = 0 or out-degree = 0 at one end)
shorter than the mean read length. Tips are traced by following single-degree
nodes until a branch point is reached. A traced tip is removed only when both
conditions hold:

- the traced chain is shorter than the mean read length
- its mean edge support is below a coverage-aware threshold

That threshold is:

$$
\max\left(0.05 \cdot \text{local avg},\; 0.25 \cdot \text{mean backbone node depth in the region}\right)
$$

where `local_avg` is the mean edge weight in a local ±500 bp window around the
tip's stored coordinates. Haplotype edges do not protect a tip from removal;
they are ignored during tip tracing and only affect downstream phasing.

### Low-weight edge pruning

For each edge, the local average weight is computed over edges within a ±500 bp
window around the stored coordinates of both endpoint nodes. Backbone nodes
contribute their exact `ref_pos`; non-backbone nodes contribute every implied
coordinate in `ref_positions`.

Backbone-to-backbone edges are removed when:

$$
   ext{edge weight} < 0.05 \cdot \text{local avg}
$$

Edges that involve at least one non-backbone node use the same regional floor
as tip pruning and are removed when:

$$
   ext{edge weight} < \max\left(0.05 \cdot \text{local avg},\; 0.25 \cdot \text{mean backbone node depth in the region}\right)
$$

This keeps pruning aggressive on weak read-derived branches without severing
weak-but-legitimate backbone continuity edges.

### Bubble popping

Bubble popping is currently disabled. The code still reports this in the
cleaning summary so debug logs make it explicit that no bubble-removal pass ran.

After each operation, the graph's adjacency lists are rebuilt to reflect
removals.

## 6. Unitig compaction

Maximal non-branching paths are collapsed into **unitigs**. In practice,
compaction starts from non-internal endpoints, extends across unique
`1-in/1-out` continuations, and then merges provisional fragments when the
cleaned DBG still implies a unique continuation across their boundary. This
keeps compaction driven by graph topology rather than node insertion order.

Singleton unitigs whose sequence length is exactly $k$ are dropped when they
have no ordinary incoming or outgoing unitig edges. Haplotype-edge-only
connectivity does not keep such a singleton in `unitig.gfa`, because haplotype
edges are phasing constraints rather than traversal edges.

Each retained unitig stores:

- The list of underlying node IDs.
- A consensus sequence (first node's k-mer plus the last character of each
  subsequent k-mer).
- The mean depth (average node depth across the chain).
- The leftmost stored coordinate (`ref_pos`) and the sorted union of all
   source-node coordinates (`ref_positions`).

Edges between unitigs are created wherever the source graph connects the last
node of one unitig to the first node of another. Edge weights are summed.

Haplotype edges are lifted from node level to unitig level. Weak haplotype
edges (weight < 5% of the minimum depth of the two connected unitigs) are
pruned, but these phasing edges are not represented as `L` lines in
`unitig.gfa`.

A DFS-based cycle detection check is run on the unitig graph. Cycles indicate
unresolvable repeat structures; if detected, the assembly aborts for that
region.

## 7. SV calling mode

SV mode reuses the cleaned unitig DAG instead of the ILP path model.

### Activation and outputs

By default, Sharda emits both haplotype and SV outputs. `--sv-only` switches
to SV-only execution, while `--hap-only` switches to haplotype-only
execution. When SV output is enabled, Sharda:

- preserves backbone-backbone edges during graph cleaning
- still emits the standard single-region graph views
- emits SV-oriented unitig views as `unitig.sv.gfa` and `unitig.sv.json`
- writes `<prefix>.sv.vcf`

With `--sv-only`, Sharda also skips haplotype flow decomposition and does not
write `<prefix>.haplotypes.fa`.

With `--hap-only`, Sharda skips SV calling, does not write `<prefix>.sv.vcf`,
and does not emit `unitig.sv.gfa` or `unitig.sv.json`.

If `--unitig-only` is also present, Sharda stops after unitig graph
construction and does not call SVs.

### Calling heuristic

The current SV caller is intentionally conservative and limited to simple
insertions and deletions.

For each backbone-only unitig with out-degree greater than one, the caller:

1. explores non-backbone alternate traversals leaving that source
2. stops each traversal at the first downstream backbone unitig where the path
   rejoins the backbone
3. enumerates all source-to-sink paths within that minimal canonical interval
4. requires exactly one all-backbone path in the interval to serve as the
   reference path
5. reconstructs reference and alternate sequences from unitig sequences using
   the stored `k-1` overlap between adjacent unitigs
6. trims shared prefix and suffix sequence
7. emits a call only if the remaining difference reduces to a simple insertion
   or deletion

### Current call ranking and collapse rules

Each emitted call carries a `SUPPORT` score, currently defined as the minimum
`mean_depth` across the non-backbone unitigs on the representative alternate
traversal.

After candidate generation, the caller applies two collapse passes:

1. exact normalized-call collapse:
   calls with identical `chrom`, normalized `POS`, `END`, `REF`, `ALT`,
   `SVTYPE`, and `SVLEN` are merged, keeping the best-supported representative
2. narrow overlapping-deletion collapse:
   overlapping canonical deletions are merged only when they share chromosome,
   `SVLEN`, retained anchor base in `ALT`, and either `SRC_REF_POS` or
   `SNK_REF_POS`

This is why the VCF can contain fewer records than the number of raw alternate
traversals in the unitig graph.

### VCF fields

`<prefix>.sv.vcf` currently contains these key INFO fields:

- `SVTYPE` — current structural-variant type (`INS` or `DEL`)
- `END` — 1-based inclusive end position of the normalized reference allele
- `SVLEN` — `ALT` length minus `REF` length
- `SUPPORT` — representative-path support score
- `SRC_UID` and `SNK_UID` — source and sink backbone unitig IDs for the
  canonical interval
- `SRC_REF_POS` and `SNK_REF_POS` — 1-based reference anchor coordinates for
  that canonical interval

The normalized VCF allele can be smaller than the full canonical interval used
to derive it. The source/sink INFO fields therefore provide additional context
about where the alternate traversal departed from and rejoined the backbone.

## 8. Flow decomposition

The goal is to decompose the unitig-level graph into *ploidy* haplotype paths
whose combined flow best explains the observed edge coverage, subject to
phasing constraints from haplotype edges.

### Source and sink identification

Flow decomposition is anchored to the region boundaries. The source is the
unitig containing the exact backbone node at the local region start, and the
sink is the unitig containing the exact backbone node at the local region end
(the backbone node at position `|ref|-k`). If either boundary anchor cannot be
resolved in the compacted graph, the implementation falls back to the previous
topology-based heuristic: use the unique in-degree-0 unitig as source and the
unique out-degree-0 unitig as sink, or the first/last such candidates by ID if
multiple exist.

### Path enumeration

All source-to-sink paths are enumerated via DFS with a backtracking visited
set, capped at 1000 paths to bound runtime.

### ILP formulation

For *P* enumerated paths and *E* unitig edges:

**Variables:**

- $f_p \geq 0$ — flow on path $p$ (continuous)
- $s_e \geq 0$ — coverage slack for edge $e$ (continuous)

**Objective:**

$$\min \sum_{e} s_e$$

**Coverage constraints** (for each edge $e$ with observed coverage $c_e$):

$$\sum_{p \in \text{paths using } e} f_p + s_e \geq c_e$$
$$\sum_{p \in \text{paths using } e} f_p - s_e \leq c_e$$

These linearise the absolute deviation: $\left|\sum_p f_p \cdot a_{pe} - c_e\right| \leq s_e$.

**Haplotype constraints** (for each haplotype edge connecting unitigs $U_1$ and
$U_2$):

$$\sum_{p \text{ passing through both } U_1 \text{ and } U_2} f_p \geq 1$$

This forces at least one path to honour the phase linkage observed from read
pairs.

### Solver

The ILP is solved by [HiGHS](https://github.com/ERGO-Code/HiGHS) with a
configurable time limit (default 60 s). If no optimal solution is found, a
fallback distributes average coverage equally across the first *ploidy* paths.

### Output extraction

Paths with flow ≥ 0.5 are kept, sorted by descending flow. Each path's
sequence is the concatenation of its unitig sequences. In debug mode, the same
extracted paths can also be persisted to `flow_paths.json`, which records each
path's flow, ordered unitigs, and per-unitig DBG node membership.

## 8. Coordinate handling in parallel mode

In whole-genome mode, each target region is assembled independently:

1. Reads overlapping the region (with configurable flanking padding) are
   extracted from the coordinate-sorted BAM via htslib index queries.
2. Extracted reads are name-sorted in memory and written to a temporary BAM.
3. The reference subsequence for the region is fetched via `faidx`.
4. Tandem repeats from the genome-wide BED are filtered and converted to local
   (0-based region-relative) coordinates.
5. The standard single-region pipeline runs.
6. Output coordinates are adjusted back to genomic space.

Regions are distributed across threads via an atomic work-stealing index,
with one `std::thread` per requested thread.
