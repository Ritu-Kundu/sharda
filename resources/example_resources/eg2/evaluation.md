# Example 2 Evaluation

This example was re-evaluated on 2026-05-16 using the current `./build/sharda`
binary from the extracted region in
`resources/example_resources/eg2/example_region.fasta`.

## Configuration

- Region: `chr5:70954000-70956000`
- Sample model: heterozygous deletion with no truth small variants
- Deleted interval: `[70954500, 70955400)` relative to chr5 coordinates
- Deleted segment length: 900 bp
- Reference-like haplotype length: 2001 bp
- Deleted haplotype length: 1101 bp
- Read model: paired-end 150 bp reads
- Sequencing depth: 10x total, split as approximately 5x per haplotype
- Read simulator: `wgsim`
- Extra simulated variants: disabled (`-r 0 -R 0 -X 0`)
- Sequencing error rate: `wgsim -e 0.01`
- Aligner: `minimap2 -ax sr`
- Assembler: `./build/sharda -d -p 2 -k 45`

## Generated Inputs

- Truth reference-like haplotype:
  `resources/example_resources/eg2/example_region_ref.fa`
- Truth deleted haplotype:
  `resources/example_resources/eg2/example_region_del.fa`
- Diploid truth FASTA: `resources/example_resources/eg2/example_sample_truth.fa`
- Truth haplotypes aligned to the example reference for IGV:
  `resources/example_resources/eg2/results/example_sample_truth_vs_ref.bam`
- Simulated reads: `resources/example_resources/eg2/reads/example_reads_R1.fq`
  and `resources/example_resources/eg2/reads/example_reads_R2.fq`
- Coordinate-sorted BAM: `resources/example_resources/eg2/results/example_reads.coord.bam`
- Name-sorted BAM for Sharda:
  `resources/example_resources/eg2/results/example_reads.namesorted.bam`

## Input Validation

- Simulated read pairs: 53
- Total mapped reads in coordinate-sorted BAM: 110 / 110
- Primary mapped reads: 106 / 106
- Properly paired reads: 74 / 106 primary reads (`69.81%`)

## Sharda Run

Command:

```bash
./build/sharda \
  -d \
  -k 45 \
  -r resources/example_resources/eg2/example_region.fasta \
  -b resources/example_resources/eg2/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg2/results/sharda_eg2_k45_componentfix
```

Observed result:

- Reference loaded: 2001 bp
- Added read pairs: 53
- Graph after read addition: 6025 nodes, 6097 edges, 16 haplotype edges
- Cleaned graph: 6025 nodes, 2004 edges
- Unitig graph: 8 unitigs, 10 edges, 3 haplotype edges
- Final status: four haplotype FASTA records and two SV calls were produced

Sharda wrote graph, haplotype, and SV outputs for the current build:

- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.raw.gfa`
- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.clean.gfa`
- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.unitig.gfa`
- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.unitig.sv.gfa`
- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.unitig.sv.json`
- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.haplotypes.fa`
- `resources/example_resources/eg2/results/sharda_eg2_k45_componentfix.sv.vcf`

Observed haplotype FASTA records:

- 2177 bp
- 1321 bp
- 1321 bp
- 2242 bp

Observed SV calls in the VCF: 2 records

The verified fresh debug rerun produced
`resources/example_resources/eg2/results/sharda_eg2_k45_componentfix_debug/`,
containing the raw, cleaned, and unitig graph bundles in both JSON and
compatibility-GFA form. In that bundle, the previously discussed weak branch
around node 5153 is removed by the current cleaner.

## Evaluation Outcome

For the current `eg2` sample, the new internal alternate-component pruning pass
collapses a large fraction of the weak 1x read-only structure, reducing the
cleaned graph from 6097 to 2004 edges and the unitig graph to 8 unitigs in the
verified fresh rerun. The current build still emits a haplotype FASTA and an SV
VCF, but the assembled paths remain longer than the 2001 bp reference and are
not yet a faithful diploid reconstruction.