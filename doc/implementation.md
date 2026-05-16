# Implementation

This document describes the code architecture, module-by-module design, data
flow, and build system of Sharda.

For user-facing debug workflow and artifact format details, see
`doc/debugging.md`.

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
Node             { id, kmer, ref_pos, ref_positions, is_backbone, depth }
Edge             { from, to, weight }
HaplotypeEdge    { from_node, to_node, weight }
Unitig           { id, node_ids, sequence, mean_depth,
                   backbone_node_count, read_node_count }
HaplotypePath    { unitig_ids, sequence, flow }
StructuralVariantCall { chrom, pos, end, id, ref, alt, sv_type, sv_len,
                        filter, info_fields }
```

Additional enums model pipeline mode and unitig provenance:

- `ExecutionMode` — `Haplotype`, `Sv`, or `Both`
- `UnitigSupportClass` — `BackboneOnly`, `Mixed`, or `ReadOnly`

`ref_pos` remains the primary scalar coordinate used by existing consumers.
For backbone nodes it is the exact backbone position. For non-backbone nodes it
is the minimum element of `ref_positions`, which stores every implied local
coordinate that has reused that read node.

### `src/graph/dbg.h / dbg.cpp` — Positional de Bruijn graph

The `DBG` class stores:

- `nodes_` — flat vector of `Node` objects.
- `edges_` — flat vector of `Edge` objects.
- `fwd_adj_` / `rev_adj_` — per-node adjacency lists (edge indices).
- `pos_to_node_` — map from `ref_pos` to node ID (backbone lookup).
- `kmer_to_node_` — map from k-mer string to node ID (read-node lookup).
- `backbone_kmer_to_nodes_` — map from backbone k-mer string to all backbone
  node IDs carrying that k-mer, used for nearest-position reuse.
- `tr_to_nodes_` — map from TR ID to the set of backbone node IDs inside that TR.
- `node_removed_` / `edge_removed_` — soft-delete flags.

Key operations:

| Method | Description |
| ------ | ----------- |
| `add_backbone_node(kmer, ref_pos)` | Positional node, indexed by `ref_pos` |
| `add_read_node(kmer)` | Hash-based node, indexed by `kmer` string |
| `closest_backbone_node_for_kmer(kmer, implied_ref_pos)` | Returns the backbone occurrence of `kmer` nearest to the implied coordinate |
| `add_node_ref_pos(node_id, ref_pos)` | Records an implied coordinate on a reused node |
| `add_edge(from, to)` | Increments weight if edge exists, else creates |
| `add_haplotype_edge(from, to)` | Increments weight if exists |
| `remove_node(id)` / `remove_edge(idx)` | Soft-delete |
| `rebuild_adjacency()` | Reconstructs adjacency from non-removed edges |
| `backbone_node_at(ref_pos)` | O(1) backbone node lookup by exact position |
| `find_read_node(kmer)` | O(1) read-node lookup by k-mer |

### `src/graph/backbone.h / backbone.cpp` — Backbone builder

`build_backbone(ref_seq, k, trs, graph)`:

1. Extracts k-mers from the reference at positions 0..len-k.
2. Creates one backbone node per position.
3. Adds edges between consecutive positions.
4. Registers each backbone node with overlapping TRs.

### `src/graph/unitig_graph.h / unitig_graph.cpp` — Unitig compaction

`UnitigGraph::build(const DBG& source)`:

1. Identifies non-internal nodes (in-degree ≠ 1 or out-degree ≠ 1).
2. Traces backward and forward from each non-internal node along unique linear
  continuations so compaction is not sensitive to node ID order.
3. Handles isolated all-internal chains as a second pass.
4. Merges provisional unitigs when the source DBG still implies a unique
  predecessor/successor continuation across the unitig boundary.
5. Builds unitig-level edges and haplotype edges by mapping node-level
   connections through the `node_to_unitig_` map.
6. Prunes weak haplotype edges (< 5% of min endpoint coverage).
7. Drops singleton unitigs with no ordinary unitig-edge incidence. Haplotype
  edges alone do not keep a singleton because they are not emitted in
  `unitig.gfa`.
8. Records how many backbone-derived and read-only DBG nodes contributed to
  each unitig so downstream SV serializers and callers can distinguish
  backbone-only, mixed, and read-only paths.

`UnitigGraph::detect_cycles()` runs iterative DFS with three-colour marking
(white/gray/black).

`UnitigGraph::k()` preserves the source DBG k-mer size so downstream path
reconstruction can append unitig suffixes using a `k-1` overlap.

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

### `src/io/debug_artifacts.h`

Thin orchestration layer for persisted debug outputs.

- `ensure_debug_output_dir(config)` — creates the artifact directory and emits
  a lightweight static `viewer.html` placeholder.
- `write_dbg_debug_artifacts(config, stage_name, graph)` — writes a stage's DBG
  artifacts (`<stage>.gfa`, `<stage>.json`) plus `manifest.json`.
- `write_unitig_debug_artifacts(config, stage_name, graph)` — writes the same
  artifact set for the compacted unitig graph.
- `write_flow_path_artifacts(config, unitig_graph, paths)` — writes
  `flow_paths.json` for extracted ILP paths, including only the flow value and
  ordered unitig IDs for each path.

SV-oriented unitig artifacts are currently written by the unitig-stage callers
in `main.cpp` and `region_assembler.cpp` rather than via a separate manifest
entry.

### `src/io/gfa_writer.h / gfa_writer.cpp`

- `write_gfa(path, graph)` — GFA1 output from a `DBG`.
- `write_dbg_json(path, graph)` — structured JSON snapshot from a `DBG`.
- `write_unitig_gfa(path, unitig_graph)` — GFA1 output from a `UnitigGraph`.
- `write_unitig_json(path, unitig_graph)` — structured JSON snapshot from a
  `UnitigGraph`, including aggregated unitig reference coordinates.
- `write_sv_unitig_gfa(path, unitig_graph)` — additive SV-oriented GFA output
  with support-class tags and color hints.
- `write_sv_unitig_json(path, unitig_graph)` — additive SV-oriented JSON output
  carrying support class and provenance counts.
- `write_vcf(path, calls, source)` — writes a minimal VCF v4.3 file for the
  current set of SV calls.
- `write_flow_path_artifacts(config, paths)` — structured JSON artifact for
  ILP output paths.

For DBG artifacts, serializer-visible node IDs are emitted in a stable order so
repeated runs on the same input produce deterministic GFA and JSON node names.

SV-oriented unitig GFA segment records currently add:

- `SC` — support class (`backbone`, `mixed`, `read`)
- `BN` — number of backbone nodes contributing to the unitig
- `RN` — number of read-only nodes contributing to the unitig
- `CL` — color hint intended for downstream graph viewers

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
2. Dispatches to `add_orr_path()` or `add_irr_path()` per read.
3. If either read is evidence: adds a haplotype edge between the pair's
   paths (last node of left read → first node of right read, ordered by
   `ref_start`).

Internal helpers:

- **`add_orr_path(read, graph)`** — derives an implied local coordinate for the
  read's first k-mer from `read.ref_start`, increments one coordinate per k-mer,
  and for each read k-mer first searches all backbone nodes carrying that k-mer.
  If multiple backbone nodes match, it reuses the one closest to the implied
  coordinate. If none match, it reuses or creates a read node keyed by k-mer
  sequence, records the implied coordinate in that node's `ref_positions`, and
  links the first divergent k-mer from backbone coordinate `N-1` when `N > 0`.

- **`add_irr_path(read, tr_id, graph)`** — extracts read k-mers, calls
  `find_and_chain_anchors`, then walks k-mer positions using anchored backbone
  nodes or hash-based read nodes. IRR anchor chaining remains the preferred
  placement path. If no anchors are found, IRR falls back to the ORR implied-
  coordinate path.

### `src/assembly/graph_cleaner.h / graph_cleaner.cpp`

`clean_graph(graph, mean_read_length, options)`:

Runs up to 10 rounds of:

1. `remove_tips(graph, max_tip_len)` — traces dead-end paths (in-degree or
  out-degree = 0) by following single-degree nodes. Removes a traced tip only
  if it is shorter than `max_tip_len` and its mean tip-edge support is below
  `max(0.05 * local_avg, 0.25 * mean_backbone_depth(graph))`. Haplotype edges
  are ignored for this decision and do not protect tips from removal.

2. `remove_weak_internal_branches(graph, max_tip_len)` — traces weak
  non-backbone alternate structure that is internal to the graph rather than a
  dead-end tip. The cleaner first handles simple linear read-only branches
  between a branching source anchor and a converging sink anchor, then falls
  back to a small non-backbone component trace that can absorb internal splits.
  A candidate branch/component is removed only when all of the following hold:
  it contains at least two non-backbone nodes, every boundary-support edge seen
  by the component has weight 1, its mean boundary support is at most 1, its
  contextual threshold `max(0.05 * local_avg, 0.25 * mean_backbone_depth(graph))`
  is not above 1, and its length stays below `max_tip_len`. This makes the pass
  target single-copy internal alternate read structure without suppressing more
  supported allelic branches.

3. `prune_low_weight_edges(graph)` — for each edge, computes the mean weight of
  all edges within ±500 bp of the edge's endpoint coordinates. A node is
  considered local to the window if its `ref_pos` or any entry in its
  `ref_positions` falls in range. Backbone-backbone edges are removed when
  `weight < 0.05 * local_avg`. Edges touching at least one non-backbone node
  instead use `weight < max(0.05 * local_avg, 0.25 * mean_backbone_depth(graph))`.

  When `options.preserve_backbone_edges=true` (currently used in SV mode),
  backbone-backbone edges are exempt from low-weight pruning so unsupported
  reference structure remains available to the SV caller.

4. Bubble popping is currently skipped. The cleaner logs
  `bubble_popping_skipped=true` in each iteration summary and does not call the
  older bubble-removal heuristic.

After each operation, `graph.rebuild_adjacency()` is called. Each iteration log
also reports `internal_branch_nodes_removed` and
`internal_branch_edges_removed` alongside the tip and low-weight edge counters.
The loop exits early if a round produces no changes.

### `src/assembly/flow_decomp.h / flow_decomp.cpp`

`flow_decomposition(unitig_graph, max_paths, anchors, time_limit_sec)`:

1. Identifies source and sink unitigs from the exact start/end backbone anchors
  when those boundary nodes survive compaction; otherwise falls back to the
  topology-based source/sink heuristic.
2. Enumerates all source-to-sink paths via DFS (capped at 1000).
3. Precomputes which edges each path uses and which unitigs each path visits.
4. Constructs an LP with HiGHS:
   - Variables: `f_p` (path flows) + `s_e` (edge coverage slack).
   - Constraints: edge coverage ± slack, haplotype edge co-occurrence ≥ 1.
   - Objective: minimise total slack.
5. Extracts paths with flow ≥ 0.5, sorted by descending flow.
6. Fallback on solver failure: equal flow distribution across first
   `max_paths` paths.

### `src/assembly/sv_caller.h` and current implementation in `region_assembler.cpp`

The initial SV caller is intentionally conservative.

`call_structural_variants(unitig_graph, chrom, coord_offset)`:

1. Computes in-degree and out-degree on the unitig DAG.
2. Chooses candidate source unitigs that are backbone-only and branch to more
  than one successor.
3. Enumerates alternate traversals starting from each candidate source and
  stops each traversal at the first downstream backbone unitig where that path
  rejoins the backbone.
4. Enumerates all source-to-sink unitig paths within that minimal canonical
  interval.
5. Requires exactly one all-backbone path to serve as the reference path.
6. Treats the canonical alternate traversal as an SV candidate when it contains
  at least one non-backbone interior unitig.
7. Reconstructs reference and alternate sequences using the stored `k-1`
  overlap between adjacent unitigs.
8. Trims shared prefix and suffix sequence and emits a call only when the
  difference reduces to a simple insertion or deletion.
9. Assigns a per-call support score equal to the minimum `mean_depth` across
  the non-backbone unitigs on the representative alternate traversal.
10. Collapses exact duplicate normalized calls, keeping the best-ranked
  representative.
11. Applies one additional conservative overlap-collapse pass for deletions
  only when two canonical calls meet all of these conditions:

- share chromosome, `SVTYPE`, and `SVLEN`
- overlap or touch in event coordinates
- retain the same anchor base in `ALT`
- share either `SRC_REF_POS` or `SNK_REF_POS`

1. Reassigns output IDs after collapsing so the returned calls are emitted as
  a compact `sv1`, `sv2`, ... sequence.

The current implementation now emits one VCF record per qualifying canonical
alternate path. This removes the earlier explosion of wider transitive
source/sink windows. Equivalent canonical calls that normalize to the same
simple indel are collapsed, keeping the strongest supported representative
traversal. Distinct overlapping canonical indels are still reported
separately, except for a narrow additional collapse step for overlapping
canonical deletions that share a source or sink boundary anchor and the same
retained anchor base. The representative path's current ranking score is
carried into the VCF as `SUPPORT`, defined as the minimum mean depth across the
non-backbone unitigs on that traversal.

### `src/assembly/region_assembler.h / region_assembler.cpp`

`assemble_region(params)` → `RegionResult`:

A self-contained, thread-safe function that runs the full single-region
assembly pipeline and returns in-memory results for one region:

1. Reads reference via `read_fasta` or `read_fasta_region`.
2. Builds backbone.
3. Iterates BAM read pairs, adding each to the graph.
4. Cleans graph.
5. Writes raw/clean debug artifacts when requested.
6. Builds unitig graph; checks for cycles.
7. Writes unitig debug artifacts, plus `unitig.sv.gfa` / `unitig.sv.json`
  when SV mode is active.
8. If `stop_after_unitig_graph` is set, returns early after artifact emission.
9. In SV mode, calls simple indels from alternate unitig paths against a unique
  backbone-only reference path and stores them in `RegionResult::sv_calls`.
10. Runs flow decomposition only when haplotype mode is enabled.
11. Converts haplotype paths into `RegionResult::haplotypes`.

`assemble_region()` does not write the final merged FASTA or VCF files. Those
top-level outputs are written by `main.cpp` after single-region execution or
after collecting all per-region `RegionResult` objects in whole-genome mode.

For SV calls, `coord_offset` is threaded into `build_indel_call()` so returned
VCF coordinates are already translated back to genomic coordinates before
`main.cpp` writes them.

`RegionParams` captures all inputs: reference path, BAM path, TRs, ploidy,
k-mer size, output prefix, debug flag, debug artifact configuration, and
optional region coordinates.

### `src/util/debug_config.h`

`DebugArtifactsConfig` centralises the first-pass debugging framework options:

- whether artifact emission is enabled
- output directory for persisted artifacts
- whether to emit GFA, JSON, and the static HTML viewer

### `src/util/kmer.h / kmer.cpp`

`extract_kmers(sequence, k)` — returns a vector of all k-mers (substrings of
length k) from the input sequence.

### `src/util/log.h`

`init_logging(debug)` — configures spdlog: debug level if `debug=true`,
info level otherwise.

### `src/main.cpp` — Entry point

Parses CLI arguments, then branches into one of three modes:

1. node lookup from an existing debug artifact directory via `--debug-dir` and
  `--debug-node`
2. whole-genome parallel assembly when `-R` is provided
3. single-region assembly otherwise

In debug mode, the single-region and per-region whole-genome paths emit a
persisted artifact bundle under `<out_prefix>_debug/` containing GFA snapshots,
JSON snapshots, a manifest, and a lightweight static HTML viewer.

Default execution currently enables both haplotype and SV output paths. That
preserves backbone-backbone edges during cleaning, emits additive
`unitig.sv.gfa` / `unitig.sv.json` outputs, and writes `<out_prefix>.sv.vcf`
with `SVTYPE`, `END`, `SVLEN`, `SUPPORT`, `SRC_UID`, `SNK_UID`,
`SRC_REF_POS`, and `SNK_REF_POS` INFO fields. The source/sink fields expose
the canonical backbone interval anchors used to derive each call. `--sv-only`
keeps that SV path but skips haplotype flow decomposition, while `--hap-only`
forces the old haplotype-only path and suppresses SV artifacts.

In single-region mode without `-d`, `main.cpp` also writes:

- `<out_prefix>.raw.gfa`
- `<out_prefix>.clean.gfa`
- `<out_prefix>.unitig.gfa`
- `<out_prefix>.unitig.sv.gfa` and `<out_prefix>.unitig.sv.json` in the
  default mode and in `--sv-only`

When `--unitig-only` is used, the program stops after unitig graph
construction and debug/unitig artifact emission. In that mode it does not call
SVs and does not write `<out_prefix>.sv.vcf`, even if `--sv-only` is also set.

- **Single-region mode** (no `-R` flag): calls the pipeline directly with the
  provided files and writes final FASTA/VCF outputs from `main.cpp`.
- **Whole-genome parallel mode** (`-R` flag):
  1. Reads target regions from the BED file.
  2. Reads the genome-wide TR BED.
  3. Spawns `j` worker threads.
  4. An `std::atomic<size_t>` index is incremented by each thread to claim the
     next unprocessed region (work-stealing pattern).
  5. Each thread calls `assemble_region()` for its claimed region after
     extracting the region BAM and reference subsequence.
  6. Results are collected and merged into one output FASTA and/or one merged
     SV VCF depending on the execution mode.

In whole-genome mode, raw/clean/unitig artifacts are written only inside each
per-region debug directory when `-d` is enabled; there is no non-debug
per-region raw/clean GFA bundle on disk.

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
| TempFileTest | 3 | FASTA round-trip, FASTA writing, and BED parsing through temporary files |
| Backbone | 2 | Backbone node/edge counts, TR registration |
| DBG | 3 | Read node creation, edge weight accumulation, haplotype edges |
| GraphCleaner | 9 | Tip pruning, regional floors, weak internal alternate branch/component pruning, and SV-mode backbone-edge preservation |
| FlowDecomposition | 1 | Boundary-anchor fallback when topology-based source/sink selection is ambiguous |
| ReadClassifier | 3 | Evidence detection for soft-clips, TR overlap, and non-evidence ORR reads |
| UnitigGraph | 4 | Linear chain compaction, branch preservation, reverse-ID-order stability, and isolated singleton dropping |
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
