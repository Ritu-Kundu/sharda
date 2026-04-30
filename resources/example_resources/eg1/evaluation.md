# Example 1 Evaluation

This example was generated on 2026-04-30 from the extracted region in
`resources/example_resources/eg1/example_region.fasta`.

## Configuration

- Region: `chr5:70954000-70956000`
- Sample model: heterozygous deletion
- Deleted interval: `[70954500, 70955400)` relative to chr5 coordinates
- Deleted segment length: 900 bp
- Reference haplotype length: 2001 bp
- Deleted haplotype length: 1101 bp
- Read model: paired-end 150 bp reads
- Sequencing depth: 10x total, split as approximately 5x per haplotype
- Read simulator: `wgsim`
- Extra simulated variants: disabled (`-r 0 -R 0 -X 0`)
- Sequencing error rate: `wgsim` default `-e 0.02`
- Aligner: `minimap2 -ax sr`
- Assembler: `./build/sharda -p 2`

## Generated Inputs

- Truth deletion haplotype: `resources/example_resources/eg1/example_region_del_het.fa`
- Diploid truth FASTA: `resources/example_resources/eg1/example_sample_truth.fa`
- Simulated reads: `resources/example_resources/eg1/results/example_reads_R1.fq`
  and `resources/example_resources/eg1/results/example_reads_R2.fq`
- Coordinate-sorted BAM: `resources/example_resources/eg1/results/example_reads.coord.bam`
- Name-sorted BAM for Sharda:
  `resources/example_resources/eg1/results/example_reads.namesorted.bam`

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
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg1/results/sharda_eg1_debug
```

Observed result:

- Reference loaded: 2001 bp
- Added read pairs: 53
- Graph after read addition: 4778 nodes, 4686 edges, 16 haplotype edges
- Graph after cleaning: 4778 nodes, 1880 edges
- Unitig graph: 2898 unitigs, 0 edges, 0 haplotype edges
- Flow decomposition: enumerated 0 source-to-sink paths
- Final status: no haplotype FASTA was produced

Sharda wrote graph debug outputs but no assembled haplotypes:

- `resources/example_resources/eg1/results/sharda_eg1.raw.gfa`
- `resources/example_resources/eg1/results/sharda_eg1.clean.gfa`
- `resources/example_resources/eg1/results/sharda_eg1.unitig.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_debug.raw.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_debug.clean.gfa`
- `resources/example_resources/eg1/results/sharda_eg1_debug.unitig.gfa`

## Evaluation Outcome

For the requested setup, the example data was created successfully and the
alignment inputs are valid, but Sharda did not recover any source-to-sink path
in the unitig graph and therefore did not emit assembled haplotypes. The
evaluation result for this exact 10x heterozygous deletion example is therefore
an assembly failure at the flow decomposition stage.