# Scoval pipeline
Command line and supplements for the project

## Introduction
We developed Scoval, a Single-cell sequencing COVerage and ALlele-based approach, to combine read-depth and phased LOH metrics, and identify neuronal CNVs that were supported by both orthogonal measurements.


## Required third party tools and packages
| software | version |
|----------|---------|
| Python   | 3.7+    |
| samtools | 1.9+    |
| bcftools |         |
| pysam    | 0.15.3+ |


## SOP to generate final callset
1. split the bam file into the single-cell bams by the cell barcode

The duplicates, low-quality and supplementary alignments should be filtered out from the original bam file. The example command lines are:
```
samtools view -bh -F 3840 input_sorted_bam > output_filtered_bam
samtools index output_filtered_bam
```
To split the bam into single-cell bams, we first sort the original bam by the barcode tag, and then split it. Users could split the bam file by their own method.  
The command lines are:
```
samtools sort -t CB -o CB_sorted_bam input_bam
python src/bamsplit.py CB_sorted_bam 
```

1. Run Ginkgo or other similar coverage-based single-cell CNV caller to generate an initial callset.  
Ginkgo command line:
```

```

2. Collect phased germline heterozygous SNPs from the same individual.

It would be at 5-columns table, including the SNP coordinates, reference and alternative allele and genotype. If the SNPs were called from 10X-linked reads, it should have another phase set column. An example can be found in `data/wg_het_snp.tsv`.

