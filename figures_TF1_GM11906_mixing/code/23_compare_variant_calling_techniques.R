library(BuenColors)
library(data.table)
library(dplyr)
library(stringr)
library(SummarizedExperiment)
library(Matrix)
"%ni%" <- Negate("%in%")

# Import variants
mgatk_vars <- read.table("../output/TF1_VMR_processed_variants_ndetect.tsv", header = FALSE, stringsAsFactors = FALSE)[,1]
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

#----------
# Visualize
#----------

# Import mgatk object
SE <- cbind(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds"), 
            readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds"))

# Find TF1 cells
rbind(read.table("../output/data1_meta.tsv", header = TRUE), 
      read.table("../output/data6_meta.tsv", header = TRUE)) %>% filter(assign == "TF1") %>%
  filter(mean_cov > 50) -> TF_cells_df

SE <- SE[,as.character(TF_cells_df$cell_id)]

ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))

# Split the heteroplasmy into forward and reverse for making these comparison plots
cov_fw <- assays(SE)[[paste0("A_counts_fw")]] + assays(SE)[[paste0("C_counts_fw")]] + assays(SE)[[paste0("G_counts_fw")]]+ assays(SE)[[paste0("T_counts_fw")]]+ 0.001
getMutMatrix_fw  <- function(letter){
  mat <- (assays(SE)[[paste0(letter, "_counts_fw")]]) / cov_fw
  rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
  ref_letter <- toupper(ref_allele)
  return(mat[ref_letter != letter & ref_letter != "N",])
}

cov_rev <- assays(SE)[[paste0("A_counts_rev")]] + assays(SE)[[paste0("C_counts_rev")]] + assays(SE)[[paste0("G_counts_rev")]]+ assays(SE)[[paste0("T_counts_rev")]] + 0.001
getMutMatrix_rev  <- function(letter){
  mat <- (assays(SE)[[paste0(letter, "_counts_rev")]]) / cov_rev
  rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
  ref_letter <- toupper(ref_allele)
  return(mat[ref_letter != letter & ref_letter != "N",])
}

# Get heteroplasmy matrix for each direction
cells_forward <- rbind(getMutMatrix_fw("A"), getMutMatrix_fw("C"), getMutMatrix_fw("G"), getMutMatrix_fw("T"))
cells_reverse <- rbind(getMutMatrix_rev("A"), getMutMatrix_rev("C"), getMutMatrix_rev("G"), getMutMatrix_rev("T"))

# aggregate specific mutations based on the above analyis
vars <- c("2573G>A", "10458C>T", "12012C>T", "3549C>A", "7399C>G", "546A>C")
df <- cbind(
  reshape2::melt(data.matrix(cells_forward[vars,])) %>% setNames(c("mut", "barcode", "forward")),
  reshape2::melt(data.matrix(cells_reverse[vars,])) %>% setNames(c("mutation", "barcode2", "reverse"))
)

# Visualize
p1 <- ggplot(df, aes(x = (forward), y = (reverse))) +
  geom_point(size = 0.01) + labs(x = "Forward heteroplasmy", y = "Reverse heteroplasmy") +
  facet_wrap(~mutation) + pretty_plot(fontsize = 8) + scale_x_continuous(limits = c(0,0.6)) + scale_y_continuous(limits = c(0,0.6))
cowplot::ggsave(p1, file = "../plots/plot_twoStrands_compare_variants.pdf", width = 3, height = 2)


