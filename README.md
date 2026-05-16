# Sharda

**S**equence **HA**plotype **R**eference-guided **D**ebruijn based **A**ssembler

Sharda is a targeted haplotype assembler for short-read sequencing data. Given
aligned reads, a reference sequence, and tandem repeat annotations, it
reconstructs phased haplotype sequences using a positional de Bruijn graph and
ILP-based flow decomposition. By default it also emits SV-oriented unitig
artifacts and calls simple indels directly from the unitig graph; `--sv-only`
keeps just the SV outputs and `--hap-only` restores the old haplotype-only
behavior.

The tandem repeat BED is optional. When omitted, Sharda runs the same assembly
pipeline without tandem repeat-aware anchor chaining.

## Features

- **Positional de Bruijn graph** — alignment-guided graph construction that
  preserves genomic coordinates, avoiding the ambiguity of hash-only DBGs.
- **Tandem repeat awareness** — anchor chaining within repeat regions for
  accurate repeat-spanning read placement.
- **Read-pair phasing** — haplotype edges link read-pair evidence across the
  graph for phase-aware path extraction.
- **ILP flow decomposition** — decomposes unitig graph coverage into haplotype
  paths with phasing constraints via the HiGHS solver.
- **Default combined output** — emits haplotypes plus SV-oriented unitig
  artifacts and simple indel calls from the same cleaned unitig graph.
- **SV-only mode** — skips haplotype flow decomposition and keeps only the SV
  outputs.
- **Haplotype-only mode** — skips SV calling and SV-oriented artifacts while
  keeping the original haplotype output path.
- **Whole-genome parallel mode** — processes many target regions concurrently
  from a coordinate-sorted BAM.

## Requirements

- C++17 compiler (GCC 8+, Clang 7+, Apple Clang 11+)
- CMake ≥ 3.14
- [htslib](https://github.com/samtools/htslib) (system-installed)
- Internet access for first build (CMake FetchContent downloads spdlog, HiGHS,
  GoogleTest)

### macOS (Homebrew)

```bash
brew install cmake htslib
```

### Ubuntu/Debian

```bash
sudo apt install cmake libhts-dev
```

## Building

```bash
git clone <repo-url> sharda
cd sharda
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

If htslib is installed in a non-standard location:

```bash
cmake -B build -DHTSLIB_DIR=/path/to/htslib/prefix
```

To skip tests:

```bash
cmake -B build -DSHARDA_BUILD_TESTS=OFF
```

## Conda Installation

This repository includes a Bioconda-style recipe in [recipe/meta.yaml](recipe/meta.yaml).

Once the package is published on Bioconda, install it with:

```bash
mamba install -c conda-forge -c bioconda sharda
```

To build the package locally from this repository:

```bash
conda install -c conda-forge -c bioconda conda-build mamba
conda build recipe
```

The recipe is configured to use conda-provided `htslib`, `highs`, and `spdlog`
instead of fetching dependencies from the network.

## Usage

### Single-region mode

Assemble haplotypes for a single target region from a name-sorted BAM:

```bash
./build/sharda \
  -r region_ref.fa \
  -b reads.namesorted.bam \
  -p 2 \
  -o output_prefix
```

Add `-t tandem_repeats.bed` to enable tandem repeat-aware IRR anchoring.

### Example data: simulated heterozygous deletion

The repository includes a small example dataset under
`resources/example_resources/eg1` built from `chr5:70954000-70956000` of
GRCh38. The example simulates a heterozygous sample with a 900 bp deletion over
`[70954500, 70955400)`, plus three SNPs, one 1 bp deletion, and one 1 bp
insertion split across the two haplotypes. The aligned reads are provided in
the format required by single-region mode.

Run Sharda on the example with:

```bash
./build/sharda \
  -d \
  -k 45 \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg1/results/sharda_eg1_k45
```

With the current build, `-d` writes the staged debug bundle under
`resources/example_resources/eg1/results/sharda_eg1_k45_debug/`.

In the current build, the standard `eg1` run emits raw, cleaned, unitig, and
SV-unitig graph outputs plus `sharda_eg1_k45.haplotypes.fa` and
`sharda_eg1_k45.sv.vcf`. The current evaluation is recorded in
`resources/example_resources/eg1/evaluation.md`.

Relevant example files:

- `resources/example_resources/eg1/example_region.fasta` — extracted reference
  region used as the assembly backbone.
- `resources/example_resources/eg1/example_region_ref_smallvars.fa` —
  reference-like haplotype carrying two SNPs and one 1 bp deletion.
- `resources/example_resources/eg1/example_region_del_smallvars.fa` — deleted
  haplotype carrying the 900 bp deletion, one SNP, and one 1 bp insertion.
- `resources/example_resources/eg1/example_sample_truth.fa` — diploid truth
  FASTA containing both variant-bearing haplotypes.
- `resources/example_resources/eg1/results/example_reads.coord.bam` —
  coordinate-sorted BAM for inspection and downstream evaluation.
- `resources/example_resources/eg1/results/example_reads.namesorted.bam` —
  name-sorted BAM consumed by Sharda single-region mode.
- `resources/example_resources/eg1/results/example_sample_truth_vs_ref.bam` —
  coordinate-sorted indexed BAM of the truth haplotypes mapped to the example
  reference for IGV visualization.

For the exact commands used to generate the example data, see
`doc/data/eg1/README.md`. For the current evaluation result on this dataset,
see `resources/example_resources/eg1/evaluation.md`.

The repository also includes a companion dataset under
`resources/example_resources/eg2` built from the same region and the same 900 bp
heterozygous deletion, but with no truth small variants and a higher simulated
sequencing error rate (`wgsim -e 0.01`).

Run Sharda on `eg2` with:

```bash
./build/sharda \
  -d \
  -k 45 \
  -r resources/example_resources/eg2/example_region.fasta \
  -b resources/example_resources/eg2/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg2/results/sharda_eg2_k45
```

In the current build, the debug run emits raw, cleaned, unitig, and SV-unitig
graph outputs plus `sharda_eg2_k45.haplotypes.fa` and `sharda_eg2_k45.sv.vcf`,
and writes the staged debug bundle under
`resources/example_resources/eg2/results/sharda_eg2_k45_debug/`. The current
evaluation is recorded in `resources/example_resources/eg2/evaluation.md`.

The latest verified local rerun of the current cleaner is summarized from
`resources/example_resources/eg2/results/sharda_eg2_k45_componentfix*` to keep
it distinct from an older local `sharda_eg2_k45_debug/` bundle.

Relevant `eg2` files:

- `resources/example_resources/eg2/example_region.fasta` — extracted reference
  region used as the assembly backbone.
- `resources/example_resources/eg2/example_region_ref.fa` — reference-like
  haplotype with no truth SNP or indel edits.
- `resources/example_resources/eg2/example_region_del.fa` — deleted haplotype
  carrying only the 900 bp deletion.
- `resources/example_resources/eg2/example_sample_truth.fa` — diploid truth
  FASTA containing the deletion-free and deletion-carrying haplotypes.
- `resources/example_resources/eg2/results/example_reads.coord.bam` —
  coordinate-sorted BAM for inspection and downstream evaluation.
- `resources/example_resources/eg2/results/example_reads.namesorted.bam` —
  name-sorted BAM consumed by Sharda single-region mode.
- `resources/example_resources/eg2/results/example_sample_truth_vs_ref.bam` —
  coordinate-sorted indexed BAM of the truth haplotypes mapped to the example
  reference for IGV visualization.

For the exact commands used to generate `eg2`, see `doc/data/eg2/README.md`.
For the current evaluation result on this dataset, see
`resources/example_resources/eg2/evaluation.md`.

### Whole-genome parallel mode

Assemble many target regions in parallel from a coordinate-sorted, indexed BAM:

```bash
./build/sharda \
  -r reference.fa \
  -b reads.sorted.bam \
  -R target_regions.bed \
  -p 2 \
  -j 8 \
  -o output_prefix
```

Add `-t tandem_repeats.bed` if genome-wide tandem repeat annotations are
available.

This requires:

- An indexed reference FASTA (`.fai` file alongside `reference.fa`)
- A coordinate-sorted BAM with an index (`.bai` file)

### SV outputs

Default execution writes both haplotype and SV outputs:

```bash
./build/sharda \
  -r region_ref.fa \
  -b reads.namesorted.bam \
  -p 2 \
  -o output_prefix
```

This produces `<prefix>.haplotypes.fa`, `<prefix>.sv.vcf`, and the SV-oriented
unitig views alongside the standard graph outputs.

### SV-only mode

Call structural variants from the unitig graph without producing haplotype
FASTA output:

```bash
./build/sharda \
  --sv-only \
  -r region_ref.fa \
  -b reads.namesorted.bam \
  -k 45 \
  -o output_prefix
```

Whole-genome SV-only mode uses the same flag with `-R`:

```bash
./build/sharda \
  --sv-only \
  -r reference.fa \
  -b reads.sorted.bam \
  -R target_regions.bed \
  -j 8 \
  -o output_prefix
```

### Haplotype-only mode

Use `--hap-only` to keep the pre-change behavior: emit haplotypes without SV
calls or SV-oriented unitig artifacts.

```bash
./build/sharda \
  --hap-only \
  -r region_ref.fa \
  -b reads.namesorted.bam \
  -p 2 \
  -o output_prefix
```

Current SV-output behavior:

- Default mode and `--sv-only` preserve backbone-backbone edges during graph
  cleaning.
- Still emits the standard single-region graph views (`raw.gfa`, `clean.gfa`,
  `unitig.gfa`) in non-debug single-region runs.
- Adds SV-oriented unitig views as `unitig.sv.gfa` and `unitig.sv.json`.
- Writes `<prefix>.sv.vcf` with `SVTYPE`, `END`, `SVLEN`, `SUPPORT`,
  `SRC_UID`, `SNK_UID`, `SRC_REF_POS`, and `SNK_REF_POS`.
- Calls sequence-resolved simple insertions and deletions only.

Important distinction:

- Default execution emits both haplotype and SV outputs.
- `--sv-only` disables haplotype flow decomposition and emits SV outputs only.
- `--hap-only` disables SV calling and SV-oriented unitig outputs while keeping
  haplotype flow decomposition.
- `--unitig-only` stops before both haplotype decomposition and SV calling.
  If `--unitig-only` is present, `<prefix>.sv.vcf` is not written even if
  `--sv-only` is also supplied.

### Options

| Flag | Description | Default |
| ---- | ----------- | ------- |
| `-r` | Reference FASTA | required |
| `-b` | BAM file | required |
| `-t` | Tandem repeat BED | optional |
| `-p` | Ploidy | required |
| `-R` | Target regions BED (enables parallel mode) | — |
| `-j` | Number of threads (with `-R`) | 1 |
| `-f` | Flanking padding in bp (with `-R`) | 1000 |
| `-k` | k-mer size | 121 |
| `-o` | Output prefix | `sharda_out` |
| `--unitig-only` | Stop after unitig graph construction; skip ILP and haplotype FASTA output | off |
| `--sv-only` | Disable haplotype decomposition and emit SV outputs only | off |
| `--hap-only` | Disable SV calling and SV-oriented artifacts; emit haplotypes only | off |
| `-d` | Enable debug logging and GFA output | off |
| `--trace-read` | Persist a trace for a specific read name in debug mode | repeatable |
| `--trace-locus` | Persist a trace for a local reference interval in debug mode | repeatable |
| `--debug-dir` | Inspect an existing debug artifact directory | — |
| `--debug-node` | Look up a node or unitig segment by name | — |
| `--debug-stage` | Restrict lookup to `raw`, `clean`, or `unitig` | all stages |

### Output

- `<prefix>.haplotypes.fa` — assembled haplotype sequences.
  In whole-genome mode, contig names embed region coordinates
  (e.g., `chr1:10000-20000_hap1_flow30`).
- `<prefix>.sv.vcf` — simple indel calls from the unitig graph.
  Written in the default mode and in `--sv-only`, unless `--unitig-only` or
  `--hap-only` is set.
- In single-region non-debug runs: `<prefix>.raw.gfa`, `<prefix>.clean.gfa`,
  and `<prefix>.unitig.gfa`.
- In the default mode and in `--sv-only`: `<prefix>.unitig.sv.gfa` and
  `<prefix>.unitig.sv.json`.
- With `-d`: a debug artifact directory named `<prefix>_debug/`.
  In single-region mode it contains `raw.gfa`, `clean.gfa`, `unitig.gfa`,
  `unitig.sv.gfa` and `unitig.sv.json` in the default mode and in `--sv-only`, `raw.json`,
  `clean.json`, `unitig.json`, `manifest.json`, `viewer.html`,
  `flow_paths.json` after ILP path extraction, and optionally
  `read_traces.json` and `locus_traces.json` when `--trace-read` or
  `--trace-locus` are used. `unitig.gfa` segment lines carry aggregated
  coordinate tags (`RP`, `RPS`), and `unitig.json` includes per-unitig
  `ref_pos` and `ref_positions` fields. `flow_paths.json` is not written when
  `--unitig-only` is used or when no ILP paths are extracted.
- In parallel mode, per-region debug output is written to
  `<prefix>_debug/<region>/`.

### Debugging

For a fuller debugging workflow and artifact format description, see
[doc/debugging.md](doc/debugging.md).

Generate debug artifacts for a single-region run:

```bash
./build/sharda \
  -d \
  -k 55 \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o /tmp/sharda_debug_demo
```

This writes `/tmp/sharda_debug_demo_debug/` with both compatibility GFA files
and structured JSON snapshots for the raw, cleaned, and unitig graphs, plus
`flow_paths.json` after ILP extraction. The unitig artifacts preserve the
leftmost coordinate represented by each unitig and the full sorted set of
member-node coordinates.

Inspect the extracted ILP paths and their unitig composition from an existing
debug bundle:

```bash
python -m json.tool /tmp/sharda_debug_demo_debug/flow_paths.json
```

Each path entry includes the estimated `flow` and the ordered `unitig_ids` for
that extracted ILP path.

Look up a node by segment name from an existing debug artifact directory:

```bash
./build/sharda \
  --debug-dir /tmp/sharda_debug_demo_debug \
  --debug-node 0 \
  --debug-stage raw
```

`--debug-stage` is optional. If omitted, the command searches `raw`, then
`clean`, then `unitig`.

The current lookup path reads GFA segment (`S`) lines, so node names are the
segment IDs used in `raw.gfa`, `clean.gfa`, or `unitig.gfa`.

Trace a specific read through the raw nodes it touched and the later graph
stages:

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

This writes `/tmp/sharda_trace_demo_debug/read_traces.json`. Each traced read
record includes the raw node list, whether each node was newly created or
reused, the node's primary local coordinate plus any stored `ref_positions`,
whether it survived graph cleaning, and which unitig it maps to after
compaction.

Trace a local reference interval through the same stages:

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

`--trace-locus START:LENGTH` uses local coordinates relative to the supplied
reference sequence in single-region mode. The run writes
`/tmp/sharda_locus_demo_debug/locus_traces.json`, including the reference
subsequence for that interval, nearby raw graph nodes, whether those nodes were
removed during cleaning, their `ref_positions` when a node is reused across
multiple implied coordinates, and their unitig mapping after compaction.

## Testing

```bash
cmake --build build --target sharda_tests
cd build && ctest --output-on-failure
```

## Project structure

```text
src/
  main.cpp                    CLI entry point and pipeline orchestration
  graph/
    types.h                   Core data types (Node, Edge, AlignedRead, etc.)
    dbg.h / dbg.cpp           Positional de Bruijn graph
    backbone.h / backbone.cpp Backbone graph from reference sequence
    unitig_graph.h / .cpp     Unitig compaction and cycle detection
  io/
    fasta_reader.h / .cpp     FASTA reading (whole-file and indexed region)
    fasta_writer.h / .cpp     FASTA output
    bam_reader.h / .cpp       BAM iteration, pairing, and region extraction
    bed_reader.h / .cpp       BED parsing for TRs and target regions
    debug_artifacts.h         Debug artifact emission interface
    gfa_writer.h / .cpp       GFA1 and JSON debug graph output
  assembly/
    read_classifier.h / .cpp  ORR/IRR classification and evidence detection
    anchor_chain.h / .cpp     Anchor finding and chaining for IRRs
    read_adder.h / .cpp       Read path integration into the graph
    graph_cleaner.h / .cpp    Tip removal, low-weight pruning, bubble popping
    flow_decomp.h / .cpp      ILP-based flow decomposition via HiGHS
    region_assembler.h / .cpp Thread-safe per-region assembly pipeline
  util/
    kmer.h / kmer.cpp         k-mer extraction
    log.h                     spdlog initialisation
    debug_config.h            Debug artifact configuration
tests/
  test_all.cpp                Unit tests (GoogleTest)
doc/
  debugging.md               Debug workflow and artifact format guide
  method.md                   Algorithm description
  implementation.md           Implementation details
```

## Dependencies

| Library | Version | Purpose | Bundling |
| ------- | ------- | ------- | -------- |
| [htslib](https://github.com/samtools/htslib) | ≥ 1.10 | BAM/FASTA I/O | System |
| [spdlog](https://github.com/gabime/spdlog) | 1.15.3 | Logging | FetchContent |
| [HiGHS](https://github.com/ERGO-Code/HiGHS) | 1.9.0 | LP/ILP solver | FetchContent |
| [GoogleTest](https://github.com/google/googletest) | 1.15.2 | Unit tests | FetchContent |

## Known limitations

- Only forward-strand k-mers are used; inversions are not detected.
- Cycles in the unitig graph cause the assembly to abort for that region.
- The ILP path enumeration has a hard cap of 1000 source-to-sink paths.
- Mean read length is hard-coded at 150 bp for tip removal thresholds.

## License

MIT License — see [LICENSE](LICENSE) file for details.
