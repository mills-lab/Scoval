# Scoval pipeline
Command line and supplements for the project

## Introduction
We developed Scoval, a Single-cell sequencing COVerage and ALlele-based approach, to combine read-depth and phased LOH metrics, and identify neuronal CNVs that were supported by both orthogonal measurements.


## Required third party tools and packages
| software | version |
|----------|---------|
| Python   | 3.7+    |
| samtools | 1.9+    |
| bcftools | 1.12+   |
| pysam    | 0.15.3+ |
| bedtools | 2.29.2+ |


## SOP to generate final callset
1. split the bam file into the single-cell bams by the cell barcode.  
The duplicates, low-quality and supplementary alignments should be filtered out from the original bam file. The example command lines are:
```
samtools view -bh -F 3840 input_sorted_bam > output_filtered_bam
samtools index output_filtered_bam
```
To split the bam into single-cell bams, we first sort the original bam by the barcode tag, and then split it. Users could split the bam file by their own method.  
The command lines are:
```
samtools sort -t CB -o CB_sorted_bam input_bam
python src/bamsplit.py -i CB_sorted_bam -b barcode_list -t N_thread -o output_dir
```
Here in the barcode list file, each line is a unique barcode(CB tag) in the bam file. N_thread is the number of threads for multi-threading.

2. Run Ginkgo or other similar coverage-based single-cell CNV caller to generate an initial callset.  
We used an adaption version of Ginkgo to call CNV. Please refer to this pipeline:
https://github.com/kunalkathuria/adaptedGinkgo


In order to provide an example result of Ginkgo for demo, there is synthetic CNV call set file `demo_results/ginkgo/CNV1`. Please note that it can not generate from the `demo_data` as Ginkgo requires more squences and our demo data is just for a small genomic region.

3. Collect phased germline heterozygous SNPs from the same individual. Please use the bam files in the `demo_data` folder as an example input. There are 20 example bam files with the sequencing reads in the region of `chr1:70719280-95261033`. One cell (`GCAAACTTCTGGTCCT`) has a heterozygous deletion in the region of `chr1:76719280-89261033`. Other cells do not have CNV. There is also an example vcf file in the `demo_data` folder.
Extract the SNP information from VCF file:
```
python src/phased_hetSNPs.py -v VCF_FILE --outdir OUT_DIR 
```
It will generate two files. One is a pickle file (phased_het_snp.pkl) containing a 6-columns table to have the detailed information about the phased heterozygous SNPs (chromosome, position, reference allele, alternative allele, genotype, phase set). Another is a tsv file (phased_het_snp_pos.tsv) containing the first two columns for the next step analysis.

To run the demo data:
```
python src/phased_hetSNPs.py -v demo_data/chr1.70719280.95261033.vcf.gz --outdir demo_results/snps/
```

4. Count the number of informative reads (reads that overlap with phased het-SNP), make the 100 het-SNPs windows, and calculate the log2 ratio for each window.
The command lines are:
```
samtools mpileup -q 13 -Q 13 -l phased_het_snp_pos.tsv single_cell_bam_file > mpileup_result
python src/count_mpileup_and_make_windows.py -s phased_het_snp.pkl -m mpileup_result -w 100 -p 100 --out window_info.pkl
```

To run the demo data:
```
for file in demo_data/*.bam
do
  barcode=$(basename "$file" .bam)
  samtools mpileup -q 13 -Q 13 -l demo_results/snps/phased_het_snp_pos.tsv "$file" > demo_results/mpileup_results/"$barcode"".mpileup" 
  python src/count_mpileup_and_make_windows.py -s demo_results/snps/phased_het_snp.pkl -m demo_results/mpileup_results/"$barcode"".mpileup" -w 100 -p 100 --out demo_results/windows/"$barcode"".window.info.pkl"
done

```


5. Merge the window information from all the cells, split the tables by chromosomes, and mask windows that has small number informative reads with NA.  
The command line is:
```
python src/merge_and_split_chr.py -i input_window_info_dir -b barcode_list -a outdir_allele_count -s outdir_sum_allele_count -l outdir_log2_ratio -o outdir_abs_log2_ratio
```
To run the demo data:
```
ls demo_data/*bam|cut -f2 -d"/" | cut -f1 -d"." > demo_results/barcode_list
python src/merge_and_split_chr.py -i demo_results/windows/ -b demo_results/barcode_list -a demo_results/allele_count -s demo_results/sum_allele_count -l demo_results/log2_ratio -o demo_results/abs_log2_ratio
```


6. Generate non-CNV permutation data. We first generate 100 random shuffling CNV callsets on the whole autosome genome.  
The command line is:
```
for i in {0..100}
do
bedtools shuffle -i ginkgo_autosome_calls.bed -g hg19.autosomes.chrom.sizes > perm_dir/perm{i}.bed
done
```
Then exclude the regions overlapped with Ginkgo CNV.  
```
python src/exclude_CNV.py ginkgo_autosome_calls.bed perm_dir barcode_list all_non_CNV_perm.bed
```

For demo, we provide an example non CNV permutation file `demo_results/nonCNV.bed`. Users should generate it based on their own data.

7. Calculate the median log2 ratio across the windows within CNV regions and permutated non-CNV regions, and get the emprical p-value across all the cells.    
The command lines are:
```
# for CNV calls
python src/cal_empiricalPvalue_across_allcells_withNA.py ginkgo_autosome_calls.bed abs_log2_ratio_dir out_CNV_pkl
# for non-CNV calls
# add the header line first. The headers are "chrom", "start", "end", "bam", "CN", "barcode"
python src/cal_empiricalPvalue_across_allcells_withNA.py all_non_CNV_perm.bed abs_log2_ratio_dir out_nonCNV_pkl
```

To run the demo data:
```
header1="chrom\tstart\tend\tbam\tCN"
sed "1i$header1" demo_results/ginkgo/CNV1 > demo_results/ginkgo/ginkgo_call.bed
header2="chrom\tstart\tend\tCN\tsize\tbarcode"
python src/cal_empiricalPvalue_across_allcells_withNA.py demo_results/ginkgo/ginkgo_call.bed demo_results/abs_log2_ratio demo_results/compare/CNV.pkl

sed "1i$header2" demo_results/nonCNV.bed > demo_results/nonCNV_with_header.bed
```
To be noticed, please not use the demo data to run this section for non-CNV file as the demo data does not include these non-CNV regions.


8. Build the Gausssian mixture model (GMM) to calculate the posteror probability for each CNV call. Users could go through the jupyter notebook `src/GMM.ipynb` to generate the filtered calls. Users should change the input and output file names in the notebook.  

9. Filter the CNVs by the fraction of informative windows. Users could go through the jupyter notebook `src/filter_by_windowFrac.ipynb` to generate the filtered calls.
