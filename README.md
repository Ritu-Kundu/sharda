# Sharda

**S**equence **HA**plotype **R**eference-guided **D**ebruijn based **A**ssembler

Sharda is a targeted haplotype assembler for short-read sequencing data. Given
aligned reads, a reference sequence, and tandem repeat annotations, it
reconstructs phased haplotype sequences using a positional de Bruijn graph and
ILP-based flow decomposition.

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
GRCh38. The example simulates a heterozygous 900 bp deletion over the interval
`[70954500, 70955400)` and provides the aligned reads in the format required by
single-region mode.

Run Sharda on the example with:

```bash
./build/sharda \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg1/results/sharda_eg1
```

To keep the graph outputs and verbose logging for debugging:

```bash
./build/sharda \
  -d \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg1/results/sharda_eg1_debug
```

Relevant example files:

- `resources/example_resources/eg1/example_region.fasta` — extracted reference
  region used as the assembly backbone.
- `resources/example_resources/eg1/example_region_del_het.fa` — deleted
  haplotype truth sequence.
- `resources/example_resources/eg1/example_sample_truth.fa` — diploid truth
  FASTA containing the reference and deleted haplotypes.
- `resources/example_resources/eg1/results/example_reads.coord.bam` —
  coordinate-sorted BAM for inspection and downstream evaluation.
- `resources/example_resources/eg1/results/example_reads.namesorted.bam` —
  name-sorted BAM consumed by Sharda single-region mode.

For the exact commands used to generate the example data, see
`doc/data/eg1/README.md`. For the current evaluation result on this dataset,
see `resources/example_resources/eg1/evaluation.md`.

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
| `-d` | Enable debug logging and GFA output | off |

### Output

- `<prefix>.haplotypes.fa` — assembled haplotype sequences.
  In whole-genome mode, contig names embed region coordinates
  (e.g., `chr1:10000-20000_hap1_flow30`).
- With `-d`: GFA files for the raw, cleaned, and unitig graphs. In parallel
  mode, per-region debug output is written to `<prefix>_debug/<region>/`.

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
    gfa_writer.h / .cpp       GFA1 debug output
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
tests/
  test_all.cpp                Unit tests (GoogleTest)
doc/
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
