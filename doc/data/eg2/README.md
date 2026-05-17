# Example 2 Data Creation

This document records how the files under
`resources/example_resources/eg2` were generated for a deletion-only version of
the existing toy example.

## Inputs

- Extracted reference region:
  `resources/example_resources/eg2/example_region.fasta`
- Source reference used earlier for extraction:
  `resources/raw_resources/human_GRCh38_no_alt_analysis_set.fasta`

The extracted region is `chr5:70954000-70956000`.

## Truth Haplotype Construction

This example keeps the same heterozygous 900 bp deletion as `eg1`, but removes
all truth SNP and indel edits from both haplotypes.

- Reference-like haplotype: identical to the extracted reference region
- Deleted haplotype: `chr5:70954500-70955399 del`

The large deletion removes the half-open interval `[70954500, 70955400)`, which
is 900 bp from the 2001 bp reference region.

Generated files:

- `resources/example_resources/eg2/example_region_ref.fa`
- `resources/example_resources/eg2/example_region_del.fa`
- `resources/example_resources/eg2/example_sample_truth.fa`

## Read Simulation

Paired-end reads were generated with `wgsim` using the same read length,
insert-size settings, haplotype-specific read counts, and disabled extra
variant simulation used for `eg1`, but with a higher sequencing error rate.

Commands used:

```bash
wgsim -e 0.01 -d 500 -s 50 -N 34 -1 150 -2 150 -r 0 -R 0 -X 0 -S 11 \
  resources/example_resources/eg2/example_region_ref.fa \
  resources/example_resources/eg2/reads/ref_hap_R1.fq \
  resources/example_resources/eg2/reads/ref_hap_R2.fq

wgsim -e 0.01 -d 500 -s 50 -N 19 -1 150 -2 150 -r 0 -R 0 -X 0 -S 29 \
  resources/example_resources/eg2/example_region_del.fa \
  resources/example_resources/eg2/reads/del_hap_R1.fq \
  resources/example_resources/eg2/reads/del_hap_R2.fq
```

The final mixed FASTQs were built by prefixing read names with `ref_` or `del_`
and concatenating the two haplotype-specific read sets into:

- `resources/example_resources/eg2/reads/example_reads_R1.fq`
- `resources/example_resources/eg2/reads/example_reads_R2.fq`

## Alignment and BAM Creation

Reads were aligned back to the extracted region reference with `minimap2` and
converted into both coordinate-sorted and name-sorted BAM files.

Commands used:

```bash
minimap2 -ax sr \
  resources/example_resources/eg2/example_region.fasta \
  resources/example_resources/eg2/reads/example_reads_R1.fq \
  resources/example_resources/eg2/reads/example_reads_R2.fq | \
  samtools view -b -o resources/example_resources/eg2/results/example_reads.raw.bam -

samtools sort \
  -o resources/example_resources/eg2/results/example_reads.coord.bam \
  resources/example_resources/eg2/results/example_reads.raw.bam

samtools index resources/example_resources/eg2/results/example_reads.coord.bam

samtools sort -n \
  -o resources/example_resources/eg2/results/example_reads.namesorted.bam \
  resources/example_resources/eg2/results/example_reads.raw.bam
```

For IGV inspection of the truth haplotypes themselves, the diploid truth FASTA
was also aligned back to the extracted reference region and converted into a
coordinate-sorted indexed BAM:

```bash
minimap2 -ax asm5 \
  resources/example_resources/eg2/example_region.fasta \
  resources/example_resources/eg2/example_sample_truth.fa | \
  samtools sort -o resources/example_resources/eg2/results/example_sample_truth_vs_ref.bam -

samtools index resources/example_resources/eg2/results/example_sample_truth_vs_ref.bam
```

Generated IGV files:

- `resources/example_resources/eg2/results/example_sample_truth_vs_ref.bam`
- `resources/example_resources/eg2/results/example_sample_truth_vs_ref.bam.bai`

The `results/` directory is reserved for derived BAMs and local Sharda run
outputs. The generated FASTQ inputs live under
`resources/example_resources/eg2/reads/`, and `resources/example_resources/eg2/results/sharda_eg2_*`
is ignored in git so local assembly artifacts stay untracked.