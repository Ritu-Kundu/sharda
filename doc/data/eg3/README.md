# Example 3 Data Creation

This document records how the files under
`resources/example_resources/eg3` were generated for an overlapping-deletions
toy example.

## Inputs

- Extracted reference region:
  `resources/example_resources/eg3/example_region.fasta`
- Source reference used earlier for extraction:
  `resources/raw_resources/human_GRCh38_no_alt_analysis_set.fasta`

The extracted region is `chr5:70954000-70956000`.

## Truth Haplotype Construction

This example uses two deletion haplotypes that start at the same left
breakpoint.

- Hap 1: `chr5:70954500-70955399 del` (900 bp)
- Hap 2: `chr5:70954500-70954529 del` (30 bp)

The deleted intervals are half-open in chr5 coordinates:

- Hap 1 removes `[70954500, 70955400)`
- Hap 2 removes `[70954500, 70954530)`

Generated files:

- `resources/example_resources/eg3/example_region_hap1_del900.fa`
- `resources/example_resources/eg3/example_region_hap2_del30.fa`
- `resources/example_resources/eg3/example_sample_truth.fa`

## Read Simulation

Paired-end reads were generated with `wgsim` using 150 bp reads, the same
insert-size settings and disabled extra variant simulation used for `eg1` and
`eg2`, but with zero simulated sequencing error.

Commands used:

```bash
wgsim -e 0 -d 500 -s 50 -N 18 -1 150 -2 150 -r 0 -R 0 -X 0 \
  resources/example_resources/eg3/example_region_hap1_del900.fa \
  resources/example_resources/eg3/reads/hap1_R1.fq \
  resources/example_resources/eg3/reads/hap1_R2.fq

wgsim -e 0 -d 500 -s 50 -N 33 -1 150 -2 150 -r 0 -R 0 -X 0 \
  resources/example_resources/eg3/example_region_hap2_del30.fa \
  resources/example_resources/eg3/reads/hap2_R1.fq \
  resources/example_resources/eg3/reads/hap2_R2.fq
```

The haplotype-specific read counts were chosen to keep each haplotype near 5x
coverage despite the large length difference between the two deleted alleles.

The final mixed FASTQs were built by prefixing read names with `hap1_` or
`hap2_` and concatenating the two haplotype-specific read sets into:

- `resources/example_resources/eg3/reads/example_reads_R1.fq`
- `resources/example_resources/eg3/reads/example_reads_R2.fq`

## Alignment and BAM Creation

Reads were aligned back to the extracted region reference with `minimap2` and
converted into both coordinate-sorted and name-sorted BAM files.

Commands used:

```bash
minimap2 -ax sr \
  resources/example_resources/eg3/example_region.fasta \
  resources/example_resources/eg3/reads/example_reads_R1.fq \
  resources/example_resources/eg3/reads/example_reads_R2.fq | \
  samtools view -b -o resources/example_resources/eg3/results/example_reads.raw.bam -

samtools sort \
  -o resources/example_resources/eg3/results/example_reads.coord.bam \
  resources/example_resources/eg3/results/example_reads.raw.bam

samtools index resources/example_resources/eg3/results/example_reads.coord.bam

samtools sort -n \
  -o resources/example_resources/eg3/results/example_reads.namesorted.bam \
  resources/example_resources/eg3/results/example_reads.raw.bam
```

For IGV inspection of the truth haplotypes themselves, the diploid truth FASTA
was also aligned back to the extracted reference region and converted into a
coordinate-sorted indexed BAM:

```bash
minimap2 -ax asm5 \
  resources/example_resources/eg3/example_region.fasta \
  resources/example_resources/eg3/example_sample_truth.fa | \
  samtools sort -o resources/example_resources/eg3/results/example_sample_truth_vs_ref.bam -

samtools index resources/example_resources/eg3/results/example_sample_truth_vs_ref.bam
```

Generated IGV files:

- `resources/example_resources/eg3/results/example_sample_truth_vs_ref.bam`
- `resources/example_resources/eg3/results/example_sample_truth_vs_ref.bam.bai`

The `results/` directory is reserved for derived BAMs and local Sharda run
outputs. The generated FASTQ inputs live under
`resources/example_resources/eg3/reads/`, and `resources/example_resources/eg3/results/sharda_eg3_*`
is ignored in git so local assembly artifacts stay untracked.