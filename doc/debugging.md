# Debugging Guide

This document describes the current debugging workflow in Sharda, the artifact
bundle written by debug mode, and how to inspect graph nodes and ILP path
outputs from persisted artifacts.

The current framework captures the raw DBG, the cleaned DBG, the compacted
unitig graph, and the extracted ILP flow paths when flow decomposition runs.

## Scope

What the current debugging framework supports:

- single-region debug artifact generation with `-d`
- per-region artifact generation in whole-genome mode
- persisted GFA snapshots for compatibility with existing graph tools
- persisted JSON snapshots for structured inspection
- a simple CLI lookup for node or unitig segment names via `--debug-dir` and
  `--debug-node`
- targeted read tracing for named reads via `--trace-read`
- targeted locus tracing for local reference intervals via `--trace-locus`
- persisted ILP path tracing via `flow_paths.json`
- post-hoc read and locus inspection from persisted JSON artifacts

What it does not yet support:

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
2. open `clean.gfa` to compare the graph after tip pruning and low-weight edge pruning
3. open `unitig.gfa` to inspect the compacted graph before flow decomposition
4. open `raw.json`, `clean.json`, or `unitig.json` if you want structured
  fields rather than GFA tags
5. open `flow_paths.json` to inspect the ILP-selected paths, their flows, and
  the ordered unitigs contributing to each extracted haplotype path

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

Those segment IDs are stable serializer IDs, not raw insertion-order node IDs,
so repeated runs on the same input produce the same GFA/JSON node names.

You can also inspect persisted JSON trace artifacts without rerunning the
pipeline:

```bash
./build/sharda \
  --debug-dir /tmp/sharda_trace_demo_debug \
  --debug-read READ_NAME
```

```bash
./build/sharda \
  --debug-dir /tmp/sharda_locus_demo_debug \
  --debug-locus 100:50
```

These post-hoc commands print the matching object from `read_traces.json` or
`locus_traces.json`.

### 4. Trace a specific read

If you already know the read name you want to follow, generate artifacts with
one or more `--trace-read` flags:

```bash
./build/sharda \
  -d \
  --trace-read READ_NAME \
  -k 55 \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o /tmp/sharda_trace_demo
```

This emits `read_traces.json` in the debug directory.

Each record corresponds to a traced read mate and captures:

- read name and mate label
- ORR or IRR classification
- evidence flag and TR assignment
- the raw node path recorded while the read was added to the DBG
- whether each node was newly created or reused
- whether each raw node was later removed by graph cleaning
- the unitig ID and unitig sequence for nodes that survived compaction

### 5. Trace a specific locus

If you want to follow a local reference interval through the pre-ILP graph
stages, generate artifacts with one or more `--trace-locus` flags:

```bash
./build/sharda \
  -d \
  --trace-locus 100:50 \
  -k 55 \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o /tmp/sharda_locus_demo
```

The `START:LENGTH` pair uses local coordinates relative to the input reference
sequence in single-region mode.

This emits `locus_traces.json` in the debug directory.

Each locus record captures:

- the requested local start and length
- the corresponding global start for whole-genome per-region runs
- the reference subsequence for that interval
- raw graph nodes whose reference positions overlap the interval
- immediate neighboring raw nodes connected by one incoming or outgoing edge
- whether each captured node survived graph cleaning
- the unitig ID and unitig sequence for surviving nodes

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
flow_paths.json
read_traces.json
locus_traces.json
```

`flow_paths.json` is written only when flow decomposition runs and returns one
or more paths. It is omitted in `--unitig-only` mode.

### `raw.gfa`

DBG snapshot immediately after:

- backbone construction
- read addition
- haplotype-edge accumulation

### `clean.gfa`

DBG snapshot after graph cleaning:

- tip removal
- low-weight edge pruning
- bubble popping skipped

### `unitig.gfa`

Compacted unitig graph built from the cleaned DBG, immediately before flow
decomposition.

This file shows only ordinary unitig graph edges as `L` lines. Haplotype edges
are preserved in `unitig.json` as phasing constraints, but they are not drawn
as GFA links.

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

### `flow_paths.json`

Structured ILP output snapshot with:

- extracted path list
- per-path flow value
- ordered unitig ID list for each path

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

### `read_traces.json`

Optional artifact written only when one or more `--trace-read` flags are
provided.

It records, for each matched read mate, the node sequence that the read walked
through during read addition and the later fate of those node IDs.

### `locus_traces.json`

Optional artifact written only when one or more `--trace-locus` flags are
provided.

It records the requested reference interval, the exact reference subsequence,
the raw graph nodes around that interval, and the later fate of those node IDs.

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
weights on link lines. If a singleton unitig appears in `unitig.gfa`, it now
has at least one ordinary unitig-edge connection elsewhere in the graph;
haplotype-edge-only singletons are filtered out before emission.

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

### Read trace JSON format

`read_traces.json` currently uses this shape:

```json
{
  "reads": [
    {
      "read_name": "READ_NAME",
      "mate": "read1",
      "read_type": "ORR",
      "is_evidence": false,
      "tr_id": -1,
      "raw_nodes": [
        {
          "node_id": 0,
          "sequence": "ACG...",
          "created": false,
          "is_backbone": true,
          "ref_pos": 0,
          "tr_id": -1,
          "removed_after_clean": false,
          "unitig_id": 0,
          "unitig_sequence": "ACGT..."
        }
      ]
    }
  ]
}
```

Interpretation:

- `created=false` means the read reused an existing node
- `created=true` means the read introduced a non-backbone read node
- `removed_after_clean=true` means the node exists in the raw path but was
  removed before unitig compaction
- `unitig_id=null` means the node did not survive into the compacted graph

### Locus trace JSON format

`locus_traces.json` currently uses this shape:

```json
{
  "loci": [
    {
      "local_start": 100,
      "length": 50,
      "global_start": 100,
      "reference_sequence": "ACGT...",
      "raw_nodes": [
        {
          "node_id": 17,
          "sequence": "ACG...",
          "is_backbone": true,
          "ref_pos": 100,
          "tr_id": -1,
          "removed_after_clean": false,
          "unitig_id": 3,
          "unitig_sequence": "ACGTA..."
        }
      ]
    }
  ]
}
```

Interpretation:

- `reference_sequence` is the requested interval from the local reference
- `raw_nodes` contains overlapping reference-positioned nodes plus one-hop graph
  neighbors for local context
- `removed_after_clean=true` means the node existed after read addition but did
  not survive graph cleaning
- `unitig_id=null` means the node did not survive into the compacted graph

### Flow path JSON format

`flow_paths.json` currently uses this shape:

```json
{
  "graph_kind": "flow_paths",
  "paths": [
    {
      "path_index": 0,
      "flow": 18.0,
      "unitig_ids": [0, 4, 7]
    }
  ]
}
```

Interpretation:

- `flow` is the ILP-estimated path flow reported for that extracted haplotype
- `unitig_ids` is the ordered compacted path used to build the haplotype
  sequence
- this file is absent when `--unitig-only` skips ILP or when no paths are
  returned

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

Cleaner debug logging now reports why a traced tip was removed or rejected,
including chain length, mean tip support, local average support, and the
applied threshold. Stage-to-stage comparison is still useful, but you can now
also inspect cleaner logs to distinguish “not a tip”, “too long”, and
“support too high” cases.

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
6. rerun with `--trace-read` for any specific read you want to follow through
  the node path and later graph stages
7. rerun with `--trace-locus` for any reference interval where you want the
  local subsequence and surrounding graph context preserved as an artifact

## Future Extensions

The artifact bundle and guide are designed to support later additions without
changing the outer workflow:

- coordinate/locus lookup
- node provenance histories
- richer manifest indexing
- a real interactive viewer that reads the JSON snapshots directly
