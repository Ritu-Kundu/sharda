# Debugging Guide

This document describes the current debugging workflow in Sharda, the artifact
bundle written by debug mode, and how to inspect graph nodes from persisted
outputs.

The current framework is intentionally scoped to the pre-ILP graph pipeline.
It captures the raw DBG, the cleaned DBG, and the compacted unitig graph.

## Scope

What the current debugging framework supports:

- single-region debug artifact generation with `-d`
- per-region artifact generation in whole-genome mode
- persisted GFA snapshots for compatibility with existing graph tools
- persisted JSON snapshots for structured inspection
- a simple CLI lookup for node or unitig segment names via `--debug-dir` and
  `--debug-node`

What it does not yet support:

- coordinate or locus lookup queries
- post-ILP path tracing
- a full interactive graph UI beyond the placeholder static HTML viewer
- node history/provenance events across every mutation step

## Typical Workflow

### 1. Generate debug artifacts

Run a single-region assembly with `-d`:

```bash
./build/sharda \
  -d \
  -k 55 \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o /tmp/sharda_debug_demo
```

This writes a debug directory at:

```text
/tmp/sharda_debug_demo_debug/
```

For whole-genome mode, the top-level debug directory contains one subdirectory
per region:

```text
<prefix>_debug/
  chr1_10000-20000/
  chr1_30000-40000/
  ...
```

### 2. Inspect the emitted files

The simplest first pass is:

1. open `raw.gfa` to inspect the initial graph after read addition
2. open `clean.gfa` to compare the graph after pruning and bubble removal
3. open `unitig.gfa` to inspect the compacted graph before flow decomposition
4. open `raw.json`, `clean.json`, or `unitig.json` if you want structured
   fields rather than GFA tags

### 3. Look up a specific node

Use the lookup mode on an existing debug directory:

```bash
./build/sharda \
  --debug-dir /tmp/sharda_debug_demo_debug \
  --debug-node 0 \
  --debug-stage raw
```

If `--debug-stage` is omitted, Sharda searches `raw`, then `clean`, then
`unitig`.

Example output:

```text
stage   raw
node    0
sequence        ACGT...
tag     DP:f:12
tag     BB:i:1
tag     RP:i:0
tag     TR:i:-1
```

This lookup currently reads GFA segment (`S`) lines, so the node name must be
the segment ID used in the relevant `.gfa` file.

## Artifact Layout

In single-region mode, `-d` produces the following files in
`<prefix>_debug/`:

```text
raw.gfa
clean.gfa
unitig.gfa
raw.json
clean.json
unitig.json
manifest.json
viewer.html
```

### `raw.gfa`

DBG snapshot immediately after:

- backbone construction
- read addition
- haplotype-edge accumulation

### `clean.gfa`

DBG snapshot after graph cleaning:

- tip removal
- low-weight edge pruning
- bubble popping

### `unitig.gfa`

Compacted unitig graph built from the cleaned DBG, immediately before flow
decomposition.

### `raw.json` and `clean.json`

Structured DBG snapshots with:

- graph kind
- `k` value
- node list
- edge list
- haplotype-edge list

### `unitig.json`

Structured unitig snapshot with:

- graph kind
- unitig list
- edge list
- haplotype-edge list

### `manifest.json`

Small index file recording:

- artifact format version
- latest emitted stage
- latest graph kind
- expected stage-to-file mapping

### `viewer.html`

Static placeholder HTML page created with each debug directory. This is not yet
an interactive graph viewer; it exists so the artifact bundle already has a
stable place for a future UI.

## Artifact Formats

### GFA format

The current lookup command reads GFA because it is already a stable,
human-readable compatibility format.

For DBG snapshots, segment records are written as:

```text
S <node_id> <kmer_sequence> DP:f:<depth> BB:i:<0|1> RP:i:<ref_pos> TR:i:<tr_id>
```

Meaning of the tags:

- `DP` — node depth
- `BB` — backbone flag
- `RP` — reference position for backbone nodes, `-1` for non-backbone nodes
- `TR` — tandem-repeat ID, `-1` if none

DBG link records are written as:

```text
L <from> + <to> + <k-1>M RC:i:<weight>
```

Unitig GFA uses the unitig sequence on segment lines and unitig-level edge
weights on link lines.

### DBG JSON format

`raw.json` and `clean.json` currently use this shape:

```json
{
  "graph_kind": "dbg",
  "k": 55,
  "nodes": [
    {
      "id": 0,
      "sequence": "ACGT...",
      "ref_pos": 0,
      "tr_id": -1,
      "is_backbone": true,
      "depth": 12,
      "removed": false
    }
  ],
  "edges": [
    {"from": 0, "to": 1, "weight": 10}
  ],
  "haplotype_edges": [
    {"from": 20, "to": 91, "weight": 3}
  ]
}
```

### Unitig JSON format

`unitig.json` currently uses this shape:

```json
{
  "graph_kind": "unitig",
  "unitigs": [
    {
      "id": 0,
      "sequence": "ACGT...",
      "mean_depth": 11.5,
      "node_ids": [0, 1, 2, 3]
    }
  ],
  "edges": [
    {"from": 0, "to": 1, "weight": 8}
  ],
  "haplotype_edges": [
    {"from": 4, "to": 9, "weight": 2}
  ]
}
```

## Current Query Model

The current node lookup is intentionally simple.

- It requires an existing debug directory.
- It reads stage GFA files rather than JSON.
- It returns the first matching segment ID from the requested stage, or from
  `raw`, `clean`, and `unitig` in that order if no stage is supplied.

This keeps the first query path lightweight while the artifact schema settles.

## Interpreting Stage Differences

A practical way to debug graph behavior is to compare the same node or local
region across the three stages.

- If a node exists in `raw.gfa` but not `clean.gfa`, it was removed by graph
  cleaning.
- If a set of raw nodes appears as a single unitig sequence in `unitig.gfa`,
  the graph was compacted across that region.
- If a path exists in the raw graph but disappears in the cleaned graph, check
  edge weights and backbone status to understand whether pruning removed it.

The current framework does not yet emit explicit “why this node was removed”
events, so stage-to-stage comparison is still the main debugging method.

## Recommended Workflow For Issues

When debugging an assembly failure or an unexpected haplotype:

1. run the relevant region with `-d`
2. inspect `raw.gfa` or `raw.json` to confirm that the expected branch or node
   is present before cleaning
3. inspect `clean.gfa` or `clean.json` to see whether cleaning removed or
   simplified it
4. inspect `unitig.gfa` or `unitig.json` to see how compaction changed the
   structure before ILP
5. use `--debug-node` on any suspicious segment IDs you find while browsing the
   GFA files

## Future Extensions

The artifact bundle and guide are designed to support later additions without
changing the outer workflow:

- coordinate/locus lookup
- node provenance histories
- richer manifest indexing
- a real interactive viewer that reads the JSON snapshots directly
