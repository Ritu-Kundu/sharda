# Example 1 Evaluation

This example was re-evaluated on 2026-05-16 using the current `./build/sharda`
binary from the extracted region in
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
- Assembler: `./build/sharda -d -p 2 -k 45`

## Generated Inputs

- Truth reference-like haplotype:
  `resources/example_resources/eg1/example_region_ref_smallvars.fa`
- Truth deleted haplotype:
  `resources/example_resources/eg1/example_region_del_smallvars.fa`
- Diploid truth FASTA: `resources/example_resources/eg1/example_sample_truth.fa`
- Truth haplotypes aligned to the example reference for IGV:
  `resources/example_resources/eg1/results/example_sample_truth_vs_ref.bam`
- Simulated reads: `resources/example_resources/eg1/reads/example_reads_R1.fq`
  and `resources/example_resources/eg1/reads/example_reads_R2.fq`
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
- Graph after read addition: 2646 nodes, 2674 edges
- Cleaned graph: 2646 nodes, 2230 edges
- Unitig graph: 21 unitigs, 27 edges, 12 haplotype edges
- Cleaning iteration 0: `internal_branch_nodes_removed=0`, `internal_branch_edges_removed=0`
- Final status: six haplotype FASTA records and three SV calls were produced

Sharda wrote graph, haplotype, and SV outputs for the current build:

- `resources/example_resources/eg1/results/sharda_eg1_k45.raw.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.clean.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.unitig.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.unitig.sv.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.unitig.sv.json`
- `resources/example_resources/eg1/results/sharda_eg1_k45.haplotypes.fa`
- `resources/example_resources/eg1/results/sharda_eg1_k45.sv.vcf`

Observed haplotype FASTA records:

- 2573 bp (`flow=3.0`)
- 1673 bp (`flow=2.0`)
- 1674 bp (`flow=2.0`)
- 2573 bp (`flow=1.0`)
- 2572 bp (`flow=1.0`)
- 1672 bp (`flow=1.0`)

Observed SV calls in the VCF: 3 records

## Evaluation Outcome

For the current `eg1` sample, Sharda completes the `k=45` run and now emits a
larger set of alternate haplotype paths plus an SV VCF. The resulting output is
still not a clean diploid reconstruction of the two truth haplotypes: the
haplotype FASTA is over-fragmented and includes paths substantially longer than
the 2001 bp reference region, but the 900 bp deletion is represented in the SV
calls. In this rerun, the new internal alternate-branch pruning path did not
activate on `eg1`; the observed graph reduction still came from tip removal and
low-weight edge pruning.
