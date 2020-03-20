# Reproducing benchmark public scRNA-seq data

In this repository subfigure, we both benchmark coverages of scRNA-seq data and evaluate mgatk in identifying putative clonal variants.

For SMART-seq2 data, we reprocessed fastq files from our previous work: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115214

### Variant calling benchmarking

```
freebayes-v1.3.1 -L D1.fbin.txt -f hg19.fasta > Donor1_freebayes.vcf
freebayes-v1.3.1 -L D2.fbin.txt -f hg19.fasta > Donor2_freebayes.vcf
bcftools mpileup -Ou -f hg19.fasta -b Donor1mito_bam_list.txt | bcftools call -mv -Ob -o Donor1_bcftools.bcf
bcftools mpileup -Ou -f hg19.fasta -b Donor2mito_bam_list.txt | bcftools call -mv -Ob -o Donor2_bcftools.bcf
```


For both sets of output files, we used `bcftools view` to stream output into this `grep | awk` command

```
grep chrM | awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' 
```

in order to focus on columns needed to evluate the called variants. These

```
Donor1_bcftools.tsv
Donor1_freebayes.tsv
Donor2_bcftools.tsv
Donor2_freebayes.tsv
```

### Benchmarking coverages

To compare the four technologies, we used the following datasets:

Technology|Description|Source
--- | --- | ---
Smart-seq2|935 CD34+ differentiated cells|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115214
10x v3|10k PBMCs|https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3
10x v2|8k PBMCs|https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k
10x 5p|8k PBMCs|https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_pbmc_5gex

Bam files were preprocessed with mgatk using default settings. 