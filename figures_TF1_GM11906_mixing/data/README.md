# Description of data files

- The GM12878 LCL ChIP-seq peak files (folder gm12878_encode) were downloaded from here:
https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/

- The `*calls.tsv` files were generated from split `.bam` files available here: 
https://github.com/caleblareau/tf1_variantcall_compare

- The `merrf_tf1.fixedwidthpeaks.bed` file was created from summits (macs2) -> proatac 
https://github.com/buenrostrolab/proatac

- The `hg19-tss.bed` file comes from UCSC via
https://github.com/buenrostrolab/tss-annotation

- The `GM11906_C1*` files are the output of `mgatk call` for the one plate of GM11906 data from
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142741

- The `*.ase.tsv.gz` files used for determinging haplotypes of the two mutations were from hisnipper executions of the cellranger-atac bam files:
https://github.com/caleblareau/hisnpper

- The `heteroplasmy_and_coverage_fov09.csv` file was generated via Zack + Tongtong from the in situ experiments
