# Example 1 Evaluation

This example was generated on 2026-04-30 from the extracted region in
`resources/example_resources/eg1/example_region.fasta`.

## Configuration

- Region: `chr5:70954000-70956000`
- Sample model: heterozygous deletion plus small variants
- Deleted interval: `[70954500, 70955400)` relative to chr5 coordinates
- Deleted segment length: 900 bp
- Reference-like haplotype variants:
  `chr5:70954120 C>G`, `chr5:70954286 G>T`, `chr5:70954451 delT`
- Deleted haplotype variants:
  `chr5:70955464 T>C`, `chr5:70955731 insA`
- Reference-like haplotype length: 2000 bp
- Deleted haplotype length: 1102 bp
- Read model: paired-end 150 bp reads
- Sequencing depth: 10x total, split as approximately 5x per haplotype
- Read simulator: `wgsim`
- Extra simulated variants: disabled (`-r 0 -R 0 -X 0`)
- Sequencing error rate: `wgsim -e 0.001`
- Aligner: `minimap2 -ax sr`
- Assembler: `./build/sharda -p 2 -k 45`

## Generated Inputs

- Truth reference-like haplotype:
  `resources/example_resources/eg1/example_region_ref_smallvars.fa`
- Truth deleted haplotype:
  `resources/example_resources/eg1/example_region_del_smallvars.fa`
- Diploid truth FASTA: `resources/example_resources/eg1/example_sample_truth.fa`
- Truth haplotypes aligned to the example reference for IGV:
  `resources/example_resources/eg1/results/example_sample_truth_vs_ref.bam`
- Simulated reads: `resources/example_resources/eg1/results/example_reads_R1.fq`
  and `resources/example_resources/eg1/results/example_reads_R2.fq`
- Coordinate-sorted BAM: `resources/example_resources/eg1/results/example_reads.coord.bam`
- Name-sorted BAM for Sharda:
  `resources/example_resources/eg1/results/example_reads.namesorted.bam`

## Input Validation

- Simulated read pairs: 53
- Total mapped reads in coordinate-sorted BAM: 110 / 110
- Primary mapped reads: 106 / 106
- Properly paired reads: 84 / 106 primary reads (`79.25%`)

## Sharda Run

Command:

```bash
./build/sharda \
  -d \
  -k 45 \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg1/results/sharda_eg1_k45
```

Observed result:

- Reference loaded: 2001 bp
- Added read pairs: 53
- Final status: one haplotype FASTA record was produced (733 bp, flow 5.0)

Sharda wrote graph outputs and one assembled haplotype:

- `resources/example_resources/eg1/results/sharda_eg1_k45.raw.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.clean.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.unitig.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.haplotypes.fa`

## Evaluation Outcome

For the current `eg1` sample, Sharda still completes the `k=45` run and emits a
single 733 bp haplotype sequence. The added SNPs and 1 bp indels change the
truth haplotypes and read set, but the observed assembly remains a partial
recovery rather than a full diploid reconstruction.
