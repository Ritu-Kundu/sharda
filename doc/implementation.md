# Implementation

This document describes the code architecture, module-by-module design, data
flow, and build system of Sharda.

---

## Build system

The project uses **CMake ≥ 3.14** with `FetchContent` for most dependencies:

| Target | Type | Contents |
| ------ | ---- | -------- |
| `sharda_lib` | Static library | All `src/**/*.cpp` except `main.cpp` |
| `sharda` | Executable | `src/main.cpp`, links `sharda_lib` + `Threads::Threads` |
| `sharda_tests` | Executable | `tests/*.cpp`, links `sharda_lib` + GoogleTest |

### Dependencies

- **htslib** — found via `find_library` / `find_path`. Provide `-DHTSLIB_DIR`
  if not in system paths.
- **spdlog v1.15.3** — header-only logging, fetched via FetchContent.
- **HiGHS v1.9.0** — LP/ILP solver. Built as a static library
  (`BUILD_SHARED_LIBS=OFF`) to avoid runtime symlink issues on macOS.
- **GoogleTest v1.15.2** — fetched only when `SHARDA_BUILD_TESTS=ON` (default).

---

## Module guide

### `src/graph/types.h` — Core data types

All major structs are defined here, with no `.cpp` file. Key types:

```text
TargetRegion     { chrom, start, end }
TandemRepeat     { id, chrom, start, end }
CigarOp          enum: MATCH, INS, DEL, SOFT_CLIP, HARD_CLIP, SKIP, PAD
CigarElement     { op, length }
AlignedRead      { name, seq, qual, ref_start, ref_end, cigar, flag }
                 Helper methods: is_reverse(), is_secondary(), is_supplementary(),
                 is_proper_pair(), is_unmapped(), mate_unmapped()
ReadPair         { read1, read2 }
ReadType         enum: ORR, IRR
ReadClassification { type, is_evidence, tr_id }
Node             { id, kmer, ref_pos, is_backbone, depth }
Edge             { from, to, weight }
HaplotypeEdge    { from_node, to_node, weight }
Unitig           { id, node_ids, sequence, mean_depth }
HaplotypePath    { unitig_ids, sequence, flow }
AlignedBase      { read_pos, ref_pos, is_insert }
```

### `src/graph/dbg.h / dbg.cpp` — Positional de Bruijn graph

The `DBG` class stores:

- `nodes_` — flat vector of `Node` objects.
- `edges_` — flat vector of `Edge` objects.
- `fwd_adj_` / `rev_adj_` — per-node adjacency lists (edge indices).
- `pos_to_node_` — map from `ref_pos` to node ID (backbone lookup).
- `kmer_to_node_` — map from k-mer string to node ID (read-node lookup).
- `tr_to_nodes_` — map from TR ID to the set of backbone node IDs inside that TR.
- `node_removed_` / `edge_removed_` — soft-delete flags.

Key operations:

| Method | Description |
| ------ | ----------- |
| `add_backbone_node(kmer, ref_pos)` | Positional node, indexed by `ref_pos` |
| `add_read_node(kmer)` | Hash-based node, indexed by `kmer` string |
| `add_edge(from, to)` | Increments weight if edge exists, else creates |
| `add_haplotype_edge(from, to)` | Increments weight if exists |
| `register_tr_node(tr_id, node_id)` | Associates a backbone node with a TR |
| `remove_node(id)` / `remove_edge(idx)` | Soft-delete |
| `rebuild_adjacency()` | Reconstructs adjacency from non-removed edges |
| `node_by_pos(ref_pos)` | O(1) backbone node lookup by position |
| `node_by_kmer(kmer)` | O(1) read-node lookup by k-mer |

### `src/graph/backbone.h / backbone.cpp` — Backbone builder

`build_backbone(ref_seq, k, trs, graph)`:

1. Extracts k-mers from the reference at positions 0..len-k.
2. Creates one backbone node per position.
3. Adds edges between consecutive positions.
4. Registers each backbone node with overlapping TRs.

### `src/graph/unitig_graph.h / unitig_graph.cpp` — Unitig compaction

`UnitigGraph::build(const DBG& source)`:

1. Identifies non-internal nodes (in-degree ≠ 1 or out-degree ≠ 1).
2. Traces forward from each non-internal node along linear chains.
3. Handles isolated all-internal chains as a second pass.
4. Builds unitig-level edges and haplotype edges by mapping node-level
   connections through the `node_to_unitig_` map.
5. Prunes weak haplotype edges (< 5% of min endpoint coverage).

`UnitigGraph::detect_cycles()` runs iterative DFS with three-colour marking
(white/gray/black).

### `src/io/fasta_reader.h / fasta_reader.cpp` — FASTA I/O

- `read_fasta(path)` — reads all sequences from a FASTA file. Returns vector of
  `{name, sequence}` pairs.
- `read_fasta_region(path, chrom, start, end)` — indexed random access via
  htslib `fai_load` / `faidx_fetch_seq`. Requires a `.fai` index.

### `src/io/bam_reader.h / bam_reader.cpp` — BAM reading

- `iterate_read_pairs(bam_path, callback)` — streams through a name-sorted BAM,
  pairing reads by query name, and invokes the callback for each `ReadPair`.
  Skips secondary/supplementary/unmapped alignments. Builds `AlignedRead` from
  htslib `bam1_t`.

- `create_region_bam(bam_path, region, padding, out_path)` — extracts reads
  overlapping `chrom:start-padding..end+padding` from an indexed BAM using
  `sam_itr_querys`. Collected reads are **name-sorted in memory** via
  `std::sort` on query name, then written to the output path. This converts
  coordinate-sorted input into the name-sorted order required by
  `iterate_read_pairs`.

### `src/io/bed_reader.h / bed_reader.cpp` — BED parsing

- `read_bed(path)` — parses a BED file into `TandemRepeat` objects with
  auto-assigned IDs.
- `read_target_regions(path)` — parses a BED file into `TargetRegion` objects.
- `filter_trs_for_region(trs, region, padding)` — selects TRs overlapping the
  padded region and converts coordinates to local (0-based, region-relative)
  space: `local_start = max(0, tr.start − (region.start − padding))`.

If no TR BED is provided at the CLI, these functions are simply not called for
TR annotations and the rest of the pipeline receives an empty `trs` vector.

### `src/io/fasta_writer.h / fasta_writer.cpp`

`write_fasta(path, sequences)` — writes sequences with 80-character line
wrapping.

### `src/io/gfa_writer.h / gfa_writer.cpp`

- `write_gfa(path, graph)` — GFA1 output from a `DBG`.
- `write_unitig_gfa(path, unitig_graph)` — GFA1 output from a `UnitigGraph`.

### `src/assembly/read_classifier.h / read_classifier.cpp`

`classify_read(read, trs)` → `ReadClassification`:

1. Scans CIGAR for soft-clips ≥ 5 bp and indels ≥ 5 bp.
2. Checks SA tag (split alignment) via string search in aux data.
3. Checks improper-pair flag.
4. Tests overlap against each TR.
5. If evidence + TR overlap → IRR; else → ORR.

### `src/assembly/anchor_chain.h / anchor_chain.cpp`

`find_and_chain_anchors(read_kmers, graph, tr_id)`:

1. Builds a map from k-mer string → backbone node IDs within the given TR.
2. Filters to unique k-mers (exactly one backbone occurrence).
3. Generates candidate anchors: `(read_kmer_index, backbone_node_id)`.
4. Runs O(n²) DP with scoring: `+1` per anchor, penalty
   `|read_gap − backbone_gap|` for gap discrepancy.
5. Backtracks from the best-scoring endpoint to extract the chain.

### `src/assembly/read_adder.h / read_adder.cpp`

`add_read_pair(pair, graph, trs)`:

1. Classifies both reads.
2. Calls `cigar_walk()` to produce `AlignedBase` lists (maps read positions to
   reference positions, marking inserts).
3. Dispatches to `add_orr_path()` or `add_irr_path()` per read.
4. If either read is evidence: adds a haplotype edge between the pair's
   paths (last node of left read → first node of right read, ordered by
   `ref_start`).

Internal helpers:

- **`cigar_walk(read)`** — walks the CIGAR, emitting one `AlignedBase` per
  consumed read base. Matches emit `{read_pos, ref_pos}`, inserts emit
  `{read_pos, ref_pos, is_insert=true}`, deletes advance `ref_pos` only.

- **`add_orr_path(read, graph, aligned_bases)`** — for each aligned base, if
  `!is_insert` and a backbone node exists at `ref_pos`, uses it; otherwise
  creates a read node. Produces a path, incrementing depths and edge weights.

- **`add_irr_path(read, tr_id, graph, aligned_bases)`** — extracts read k-mers,
  calls `find_and_chain_anchors`, then walks k-mer positions using anchored
  backbone nodes or hash-based read nodes. Falls back to ORR if no anchors.

### `src/assembly/graph_cleaner.h / graph_cleaner.cpp`

`clean_graph(graph, mean_read_length)`:

Runs up to 10 rounds of:

1. `remove_tips(graph, max_tip_len)` — traces dead-end paths (in-degree or
   out-degree = 0) by following single-degree nodes. Removes paths shorter than
   `max_tip_len` (set to `mean_read_length`, default 150).

2. `prune_low_weight_edges(graph)` — for each edge, computes the mean weight of
   all edges within ±500 bp of the source node's reference position. Removes
   edges with weight < 5% of that local average.

3. `pop_bubbles(graph)` — from each node with out-degree = 2, traces both
   branches forward (up to 2k steps). If they reconverge and one branch has
   weight < 20% of the other, removes the weaker branch's non-backbone nodes.

After each operation, `graph.rebuild_adjacency()` is called. The loop exits
early if a round produces no changes.

### `src/assembly/flow_decomp.h / flow_decomp.cpp`

`flow_decomposition(unitig_graph, max_paths, time_limit_sec)`:

1. Identifies source (in-degree 0) and sink (out-degree 0) unitigs.
2. Enumerates all source-to-sink paths via DFS (capped at 1000).
3. Precomputes which edges each path uses and which unitigs each path visits.
4. Constructs an LP with HiGHS:
   - Variables: `f_p` (path flows) + `s_e` (edge coverage slack).
   - Constraints: edge coverage ± slack, haplotype edge co-occurrence ≥ 1.
   - Objective: minimise total slack.
5. Extracts paths with flow ≥ 0.5, sorted by descending flow.
6. Fallback on solver failure: equal flow distribution across first
   `max_paths` paths.

### `src/assembly/region_assembler.h / region_assembler.cpp`

`assemble_region(params)` → `RegionResult`:

A self-contained, thread-safe function that runs the full single-region
pipeline:

1. Reads reference via `read_fasta` or `read_fasta_region`.
2. Builds backbone.
3. Iterates BAM read pairs, adding each to the graph.
4. Cleans graph.
5. Builds unitig graph; checks for cycles.
6. Runs flow decomposition.
7. Writes FASTA + optional debug GFA.
8. Adjusts output coordinates by adding the region's genomic offset.

`RegionParams` captures all inputs: reference path, BAM path, TRs, ploidy,
k-mer size, output prefix, debug flag, and optional region coordinates.

### `src/util/kmer.h / kmer.cpp`

`extract_kmers(sequence, k)` — returns a vector of all k-mers (substrings of
length k) from the input sequence.

### `src/util/log.h`

`init_logging(debug)` — configures spdlog: debug level if `debug=true`,
info level otherwise.

### `src/main.cpp` — Entry point

Parses CLI arguments via `getopt`, then branches:

- **Single-region mode** (no `-R` flag): calls the pipeline directly with the
  provided files.
- **Whole-genome parallel mode** (`-R` flag):
  1. Reads target regions from the BED file.
  2. Reads the genome-wide TR BED.
  3. Spawns `j` worker threads.
  4. An `std::atomic<size_t>` index is incremented by each thread to claim the
     next unprocessed region (work-stealing pattern).
  5. Each thread calls `assemble_region()` for its claimed region after
     extracting the region BAM and reference subsequence.
  6. Results are collected and merged into a single output FASTA.

In both modes, `-t` is optional. When omitted, `main.cpp` skips `read_bed()`
and passes an empty TR list into backbone construction, read classification,
and per-region filtering.

---

## Data flow

### Single-region pipeline

```text
reference.fa ──→ read_fasta() ──→ backbone ──→ ┐
                                                ├──→ DBG
reads.bam ──→ iterate_read_pairs() ──→ classify + add_read_pair() ──→ ┘
                                                                        │
                                                        clean_graph() ←─┘
                                                              │
                                                    UnitigGraph::build()
                                                              │
                                                    flow_decomposition()
                                                              │
                                                     write_fasta() / write_gfa()
```

### Whole-genome pipeline

```text
targets.bed ──→ read_target_regions()
                       │
                       ▼
              ┌── for each region (parallel) ─────────────────────┐
              │                                                    │
              │  create_region_bam() ──→ temp name-sorted BAM      │
              │  read_fasta_region() ──→ region reference           │
              │  filter_trs_for_region() ──→ local TRs             │
              │  assemble_region() ──→ RegionResult                │
              │                                                    │
              └────────────────────────────────────────────────────┘
                       │
                       ▼
              merge results ──→ output.haplotypes.fa
```

---

## Coordinate systems

Two coordinate spaces are used:

| Space | Origin | Used in |
| ----- | ------ | ------- |
| **Genomic** | Chromosome position (0-based) | BAM records, BED files, output contig names |
| **Local** | 0-based offset within the extracted region (including padding) | Backbone `ref_pos`, internal graph operations |

In whole-genome mode, the region assembler receives local-coordinate data and
adds back the genomic offset in the output haplotype names.

TR filtering converts genomic TR coordinates to local coordinates:

```text
local_start = max(0, tr.start − (region.start − padding))
local_end   = min(region_length, tr.end − (region.start − padding))
```

---

## Threading model

Whole-genome mode uses a simple work-stealing pattern:

```cpp
std::atomic<size_t> next_region{0};
std::vector<std::thread> workers(num_threads);
for (auto& w : workers) {
    w = std::thread([&] {
        while (true) {
            size_t idx = next_region.fetch_add(1);
            if (idx >= regions.size()) break;
            results[idx] = assemble_region(params_for(idx));
        }
    });
}
for (auto& w : workers) w.join();
```

Each region runs in complete isolation — no shared mutable state beyond the
atomic index. Temporary BAM files use region-specific paths to avoid conflicts.

---

## Testing

Tests are in `tests/test_all.cpp` using GoogleTest. Current suites:

| Suite | Tests | What is covered |
| ----- | ----- | --------------- |
| Kmer | 3 | k-mer extraction: basic, short input, empty |
| TempFileTest | 1 | Temporary file creation and cleanup |
| Backbone | 2 | Backbone node/edge counts, TR registration |
| DBG | 3 | Read node creation, edge weight accumulation, haplotype edges |
| ReadClassifier | 5 | Evidence detection for soft-clips, indels, SA tags, pairs, non-evidence |
| UnitigGraph | 3 | Linear chain compaction, branching, unitig-level edge creation |
| TRFilter | 2 | Region filtering and local coordinate conversion |

Run with:

```bash
cd build && ctest --output-on-failure
```

---

## Error handling and logging

- All I/O errors (file not found, BAM open failure, FASTA index missing) are
  reported via `spdlog::error()` and cause the function to return an empty
  result.
- Per-region assembly failures in parallel mode are logged but do not abort
  other regions.
- HiGHS solver failures trigger a fallback (equal flow distribution) rather
  than aborting.
- Debug logging (`-d` flag) emits detailed per-step information: read counts,
  graph sizes after cleaning, path enumeration counts, and solver status.
