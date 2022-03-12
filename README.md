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

2. Run Ginkgo or other similar coverage-based single-cell CNV caller to generate an initial callset  
Ginkgo command line:
```

```

3. Collect phased germline heterozygous SNPs from the same individual

Extract the SNP information from VCF file:
```
python src/phased_hetSNPs.py -v VCF_FILE
```
It will generate two files. One is a pickle file (phased_het_snp.pkl) containing a 6-columns table to have the detailed information about the phased heterozygous SNPs (chromosome, position, reference allele, alternative allele, genotype, phase set). Another is a tsv file (phased_het_snp_pos.tsv) contaiing the first two columns for the next step analysis.

4. Count the number of informative reads (reads that overlap with phased het-SNP), and make the 100 het-SNPs windows  
The command lines are:
```
samtools mpileup -q 13 -Q 13 -l phased_het_snp_pos.tsv single_cell_bam_file > mpileup_result
python count_mpileup_and_make_windows.py -s phased_het_snp.pkl -m mpileup_result -w 100 -s 100 --out window_info.pkl
```

