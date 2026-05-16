# Example 3 Evaluation

This example was generated and input-validated on 2026-05-16 using the current
workspace tools and the extracted region in
`resources/example_resources/eg3/example_region.fasta`.

## Configuration

- Region: `chr5:70954000-70956000`
- Sample model: overlapping deletion haplotypes with no truth small variants
- Hap 1 deleted interval: `[70954500, 70955400)` relative to chr5 coordinates
- Hap 1 deleted segment length: 900 bp
- Hap 2 deleted interval: `[70954500, 70954530)` relative to chr5 coordinates
- Hap 2 deleted segment length: 30 bp
- Hap 1 length: 1101 bp
- Hap 2 length: 1971 bp
- Read model: paired-end 150 bp reads
- Sequencing depth: approximately 10x total, split as approximately 5x per haplotype
- Read simulator: `wgsim`
- Extra simulated variants: disabled (`-r 0 -R 0 -X 0`)
- Sequencing error rate: `wgsim -e 0`
- Aligner: `minimap2 -ax sr`
- Assembler: `./build/sharda -d -p 2 -k 45`

## Generated Inputs

- Truth haplotype with the 900 bp deletion:
  `resources/example_resources/eg3/example_region_hap1_del900.fa`
- Truth haplotype with the 30 bp deletion:
  `resources/example_resources/eg3/example_region_hap2_del30.fa`
- Diploid truth FASTA: `resources/example_resources/eg3/example_sample_truth.fa`
- Truth haplotypes aligned to the example reference for IGV:
  `resources/example_resources/eg3/results/example_sample_truth_vs_ref.bam`
- Simulated reads: `resources/example_resources/eg3/reads/example_reads_R1.fq`
  and `resources/example_resources/eg3/reads/example_reads_R2.fq`
- Coordinate-sorted BAM: `resources/example_resources/eg3/results/example_reads.coord.bam`
- Name-sorted BAM for Sharda:
  `resources/example_resources/eg3/results/example_reads.namesorted.bam`

## Input Validation

- Simulated read pairs: 51
- Total reads in coordinate-sorted BAM: 107
- Total mapped reads in coordinate-sorted BAM: 107 / 107
- Primary mapped reads: 102 / 102
- Properly paired reads: 84 / 102 primary reads (`82.35%`)

## Assembly Run

Command:

```bash
./build/sharda \
  -d \
  -k 45 \
  -r resources/example_resources/eg3/example_region.fasta \
  -b resources/example_resources/eg3/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg3/results/sharda_eg3_k45
```

Observed result:

- Final status: two haplotype FASTA records and one SV call were produced
- Unitig graph: 407 unitigs
- Flow decomposition: 2 paths extracted for ploidy 2

Sharda wrote the current outputs under the canonical prefix:

- `resources/example_resources/eg3/results/sharda_eg3_k45.haplotypes.fa`
- `resources/example_resources/eg3/results/sharda_eg3_k45.sv.vcf`
- `resources/example_resources/eg3/results/sharda_eg3_k45_debug/`

Observed haplotype FASTA records:

- 1189 bp
- 2133 bp

Observed SV calls in the VCF: 1 record

## Evaluation Outcome

The new `eg3` dataset is available as a validated overlapping-deletions input
set with indexed truth/reference FASTA files, mixed FASTQ reads, and aligned
BAMs. After rebuilding Sharda from a clean local `build/` tree on 2026-05-16,
the standard single-region `./build/sharda -d -k 45` run completes on this
sample and emits a haplotype FASTA, an SV VCF, and the staged debug bundle
under `sharda_eg3_k45_debug/`. The two assembled haplotypes still do not match
the truth lengths exactly, but `eg3` is now a verified runnable example rather
than a hanging input case.