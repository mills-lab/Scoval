# Scoval_pipeline
Command line and supplements for the project

## Introduction
We developed Scoval, a Single-cell sequencing COVerage and ALlele-based approach, to combine read-depth and phased LOH metrics, and identify neuronal CNVs that were supported by both orthogonal measurements.

## SOP to generate final callset
1. Run Ginkgo or other similar coverage-based single-cell CNV caller to generate an initial callset.  
Ginkgo command line:
```

```

2. Collect phased germline heterozygous SNPs from the same individual.

It would be at 5-columns table, including the SNP coordinates, reference and alternative allele and genotype. If the SNPs were called from 10X-linked reads, it should have another phase set column. An example can be found in `data/wg_het_snp.tsv`.

