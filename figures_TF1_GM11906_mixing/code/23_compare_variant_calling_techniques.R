library(BuenColors)
library(data.table)
library(dplyr)
library(stringr)

"%ni%" <- Negate("%in%")

mgatk_vars <- read.table("../data/VMR_processed_variants_ndetect.tsv", header = FALSE, stringsAsFactors = FALSE)[,1]
homoplasmic <- fread("../data/variant_order.txt")[[1]]

# There were 13 variants that were present in both TF1 and GM11906-- exclude these as we are comparing only 
homoplasmic_both <- c("11719G>A", "489T>C" ,  "73A>G","263A>G" ,"750A>G",   "1438A>G","2706A>G","4769A>G","8860A>G","10398A>G","15326A>G","7028C>T","14766C>T")
homoplasmic <- c(homoplasmic, homoplasmic_both)

# Filter variants
process_vcf_tsv <- function(dt, qual_val = 100){
  dt_filt <- dt %>% filter(qual >= qual_val & nchar(ref) == 1) # focus on SNVs
  # manually verified that none are multi alt alleles
  var <- paste0(as.character(dt_filt[["start"]]), as.character(dt_filt[["ref"]]), ">",  as.character(dt_filt[["alt"]]))
  var[var %ni% homoplasmic]
}

# Import variant calls
samtools_dt <- fread("../data/TF1_samtools_calls.tsv", skip = 1, header = FALSE, col.names = c("chr", "start", "ref", "alt", "qual"))
freebayes_dt <- fread("../data/TF1_freebayes_calls.tsv", skip = 1, header = FALSE, col.names = c("chr", "start", "ref", "alt", "qual"))

samtools_vars <- process_vcf_tsv(samtools_dt)
freebayes_vars <- process_vcf_tsv(freebayes_dt)

# Look at overlap between two vectors of variants
make_table <- function(vec1, vec2){
  both <- intersect(vec1, vec2)
  data.frame(both = length(both),
             vec1_only = length(vec1) - length(both), 
             vec2_only = length(vec2) - length(both))
}
make_table(mgatk_vars, samtools_vars)
make_table(mgatk_vars, freebayes_vars)

samtools_vars[samtools_vars %ni% mgatk_vars]
freebayes_vars[freebayes_vars %ni% mgatk_vars]
mgatk_vars[mgatk_vars %ni% c(samtools_vars, freebayes_vars)]


