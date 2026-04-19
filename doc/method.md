# Method

This document describes the algorithmic method used by Sharda to produce phased
haplotype sequences from aligned short reads.

The tandem repeat BED is optional. When no TR annotations are supplied, the
pipeline still runs, but all reads are effectively treated as ORR reads and the
IRR anchor-chaining path is never activated.

## Overview

Sharda follows a six-stage pipeline for each target region:

1. **Backbone construction** — build a linear chain of positional de Bruijn
   graph nodes from the reference sequence.
2. **Read addition** — thread aligned reads through the graph, choosing between
   positional (ORR) and anchor-chained (IRR) strategies.
3. **Graph cleaning** — iteratively remove tips, low-weight edges, and bubbles.
4. **Unitig compaction** — collapse maximal non-branching paths into unitigs.
5. **Flow decomposition** — solve an ILP to decompose unitig coverage into
   haplotype paths subject to phasing constraints.
6. **Output** — emit haplotype FASTA sequences (and optional GFA debug graphs).

---

## 1. Positional de Bruijn graph

A standard de Bruijn graph maps each k-mer to a single node, which is
problematic in tandem repeat regions where the same k-mer appears at multiple
genomic positions. Sharda's graph is *positional*: backbone nodes are keyed by
their reference coordinate, while non-backbone (read-derived) nodes are keyed
by k-mer string. This preserves the linear order of the reference and provides
unique anchor points inside repetitive regions.

Each node stores:

| Field | Meaning |
| ----- | ------- |
| `id` | Unique 64-bit identifier |
| `kmer` | The k-mer string (length *k*) |
| `ref_pos` | Genomic position if backbone; −1 otherwise |
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

### ORR path (alignment-guided)

For each match/mismatch operation in the CIGAR string, the corresponding
reference position is used to look up the backbone node. If the position maps
to a backbone node and the k-mer matches, the read walks along the backbone.
If there is a mismatch, a new read node is created.

Insert operations (`I`) create read nodes linked into the path. Delete
operations (`D`) skip reference positions, advancing along the backbone without
creating nodes.

This produces a path of node IDs through the graph that faithfully follows the
alignment.

### IRR path (anchor chaining)

For reads that overlap a tandem repeat, the ORR walk would place repeated
k-mers at the wrong backbone positions. Instead:

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
5. **Fallback**: if no anchors are found, fall back to the ORR path.

### Haplotype edges

After placing both reads of a pair, if either read is an evidence read, a
**haplotype edge** is added between the last node of the left read and the first
node of the right read (ordered by alignment position). Haplotype edges are not
ordinary graph edges — they do not participate in path traversal during graph
cleaning. Instead, they are used as constraints during flow decomposition to
enforce phasing.

## 5. Graph cleaning

Graph cleaning removes noise and errors via three iterative operations (up to
10 rounds, stopping when no changes are made):

### Tip removal

A **tip** is a dead-end path (in-degree = 0 or out-degree = 0 at one end)
shorter than the mean read length. Tips are traced by following single-degree
nodes until a branch point is reached. If the traced path is short, all
non-backbone nodes on it are removed.

### Low-weight edge pruning

For each edge, the local average weight is computed over edges within a ±500 bp
window around the source node's reference position. Edges with weight below 5%
of the local average are removed.

### Bubble popping

A **bubble** is detected when a node has exactly two outgoing edges whose
forward linear traces reconverge at the same node (or share a common successor)
within 2*k* steps. If one branch has weight below 20% of the other, the weaker
branch's non-backbone nodes are removed.

After each operation, the graph's adjacency lists are rebuilt to reflect
removals.

## 6. Unitig compaction

Maximal non-branching paths (chains of nodes with in-degree 1 and out-degree 1)
are collapsed into **unitigs**. Each unitig stores:

- The list of underlying node IDs.
- A consensus sequence (first node's k-mer plus the last character of each
  subsequent k-mer).
- The mean depth (average node depth across the chain).

Edges between unitigs are created wherever the source graph connects the last
node of one unitig to the first node of another. Edge weights are summed.

Haplotype edges are lifted from node level to unitig level. Weak haplotype
edges (weight < 5% of the minimum depth of the two connected unitigs) are
pruned.

A DFS-based cycle detection check is run on the unitig graph. Cycles indicate
unresolvable repeat structures; if detected, the assembly aborts for that
region.

## 7. Flow decomposition

The goal is to decompose the unitig-level graph into *ploidy* haplotype paths
whose combined flow best explains the observed edge coverage, subject to
phasing constraints from haplotype edges.

### Source and sink identification

The unique unitig with in-degree 0 is the source; the one with out-degree 0 is
the sink. If multiple candidates exist, the first source and last sink (by ID)
are selected.

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
sequence is the concatenation of its unitig sequences.

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
