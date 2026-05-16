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
  -o resources/example_resources/eg2/results/sharda_eg2_k45
```

Observed result:

- Final status: four haplotype FASTA records and two SV calls were produced
- Unitig graph: 651 unitigs
- Flow decomposition: 2 paths extracted for ploidy 2

Sharda wrote graph, haplotype, and SV outputs for the current build:

- `resources/example_resources/eg2/results/sharda_eg2_k45.haplotypes.fa`
- `resources/example_resources/eg2/results/sharda_eg2_k45.sv.vcf`
- `resources/example_resources/eg2/results/sharda_eg2_k45_debug/`

Observed haplotype FASTA records:

- 2177 bp
- 1321 bp
- 1321 bp
- 2242 bp

Observed SV calls in the VCF: 2 records

The current debug rerun produced
`resources/example_resources/eg2/results/sharda_eg2_k45_debug/`, containing the
raw, cleaned, and unitig graph bundles in both JSON and compatibility-GFA form.

## Evaluation Outcome

For the current `eg2` sample, Sharda completes the canonical `-d -k 45` run and
emits a haplotype FASTA, an SV VCF, and the staged debug bundle under
`sharda_eg2_k45_debug/`. The assembled paths remain longer than the 2001 bp
reference and are not yet a faithful diploid reconstruction, but the rerun is
now recorded against the standard `sharda_eg2_k45*` prefix rather than the
older local `componentfix` outputs.