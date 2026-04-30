# Example 1 Data Creation

This document records how the files under
`resources/example_resources/eg1` were generated for the simulated deletion
example.

## Inputs

- Extracted reference region:
  `resources/example_resources/eg1/example_region.fasta`
- Source reference used earlier for extraction:
  `resources/raw_resources/human_GRCh38_no_alt_analysis_set.fasta`

The extracted region is `chr5:70954000-70956000`.

## Deleted Haplotype Construction

The example uses a heterozygous deletion covering the half-open interval
`[70954500, 70955400)`, which removes 900 bp from the 2001 bp reference region.

The deleted haplotype FASTA was created from the extracted region by removing
bases `500..1399` in 0-based coordinates relative to the extracted sequence.

Generated files:

- `resources/example_resources/eg1/example_region_del_het.fa`
- `resources/example_resources/eg1/example_sample_truth.fa`

## Read Simulation

Paired-end reads were generated with `wgsim` using 150 bp reads, default insert
size settings, and no extra simulated SNP/indel variation beyond sequencing
errors.

Commands used:

```bash
wgsim -e 0.02 -d 500 -s 50 -N 34 -1 150 -2 150 -r 0 -R 0 -X 0 -S 11 \
  resources/example_resources/eg1/example_region.fasta \
  resources/example_resources/eg1/results/ref_hap_R1.fq \
  resources/example_resources/eg1/results/ref_hap_R2.fq

wgsim -e 0.02 -d 500 -s 50 -N 19 -1 150 -2 150 -r 0 -R 0 -X 0 -S 29 \
  resources/example_resources/eg1/example_region_del_het.fa \
  resources/example_resources/eg1/results/del_hap_R1.fq \
  resources/example_resources/eg1/results/del_hap_R2.fq
```

The final mixed FASTQs were built by prefixing read names with `ref_` or `del_`
and concatenating the two haplotype-specific read sets into:

- `resources/example_resources/eg1/results/example_reads_R1.fq`
- `resources/example_resources/eg1/results/example_reads_R2.fq`

## Alignment and BAM Creation

Reads were aligned back to the extracted region reference with `minimap2` and
converted into both coordinate-sorted and name-sorted BAM files.

Commands used:

```bash
minimap2 -ax sr \
  resources/example_resources/eg1/example_region.fasta \
  resources/example_resources/eg1/results/example_reads_R1.fq \
  resources/example_resources/eg1/results/example_reads_R2.fq | \
  samtools view -b -o resources/example_resources/eg1/results/example_reads.raw.bam -

samtools sort \
  -o resources/example_resources/eg1/results/example_reads.coord.bam \
  resources/example_resources/eg1/results/example_reads.raw.bam

samtools index resources/example_resources/eg1/results/example_reads.coord.bam

samtools sort -n \
  -o resources/example_resources/eg1/results/example_reads.namesorted.bam \
  resources/example_resources/eg1/results/example_reads.raw.bam
```

## Assembly Run

Sharda was run in single-region mode with ploidy 2:

```bash
./build/sharda \
  -d \
  -r resources/example_resources/eg1/example_region.fasta \
  -b resources/example_resources/eg1/results/example_reads.namesorted.bam \
  -p 2 \
  -o resources/example_resources/eg1/results/sharda_eg1_debug
```

This run produced graph GFA outputs but no haplotype FASTA. The detailed run
result is summarized in `resources/example_resources/eg1/evaluation.md`.