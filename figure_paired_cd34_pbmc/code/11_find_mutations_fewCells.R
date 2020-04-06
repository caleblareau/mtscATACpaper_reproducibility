library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
source("../../global_functions/variant_calling.R")
"%ni%" <- Negate("%in%")

called_mutations <- readRDS("../output/filteredCD34_mgatk_calls.rds")

# Import mgatk object
import_mgatk <- function(library){
  df <- fread(paste0("../output/barcode_qc/",library,".barcode_qc.tsv")) %>%
    filter(keep)
  se <- readRDS(paste0("../../../mtscATACpaper_large_data_files/source/mgatk_output/",library,"_v12-mtMask_mgatk.rds"))
  se_filt <- se[,as.character(df$sample)]
  colnames(se_filt) <- paste0(library, "_", colnames(se_filt))
  se_filt
}

#--------------------
# Call CD34 mutations in small numbers of cells
#--------------------
SE_CD34 <- cbind(import_mgatk("CD34_G10"), import_mgatk("CD34_H8"))
mut_se_CD34 <- call_mutations_mgatk(SE_CD34)
misc_df_CD34 <- data.frame(rowData(mut_se_CD34))

df <- misc_df_CD34 %>% filter(n_cells_conf_detected > 0 & n_cells_conf_detected <= 3 & n_cells_over_20 > 0) %>% filter(n_cells_over_5 <= 5) 

# See if these newly-identified variants are present in the existing calls
intersect(df$variant, rowData(called_mutations)$variant)
intersect(df$variant, rowData(called_mutations)$variant[rowData(called_mutations)$called_in_CD34])

df <- df[df$variant %ni% rowData(called_mutations)$variant,]
dim(df)

called_variants <- df$variant

#--- 
# Reuse code to make mutational signature
#---

# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

# Process 3 digit signature based on letters
ref_all <- fread("../../figure_invitro_CD34/data/chrM_refAllele.txt")
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)

# Annotate with called variants

ref_all_long$called <- ref_all_long$variant %in% called_variants

# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)

prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "../plots/all_mito_signature_rare_variants.pdf", width = 3, height = 1.8)


# --
# Visualize the distributions of  mean allele frequency as a sanity check
# --

af_called_in_mgatk <- rowData(called_mutations)[rowData(called_mutations)$called_in_CD34,"mean"]
af_called_rare <- df$mean

apdf <- data.frame(
  AF = c(af_called_in_mgatk, af_called_rare),
  what = c(rep("m", length(af_called_in_mgatk)), rep("r", length(af_called_rare)))
)

pDens <- ggplot(apdf, aes(x = AF*100, fill = what)) + 
  geom_density(alpha = 0.5) +
  scale_fill_manual(values =c ("firebrick", "dodgerblue3")) +
  scale_y_continuous(expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Pseudobulk % Heteroplasmy", y = "empirical density") +
  theme(legend.position = "none")
cowplot::ggsave2(pDens, file = "../plots/meanAF_density.pdf", width = 1.8, height = 1.8)

# --
# Find the top heteroplasmy per mutation for well covered (min 20x) cells
# --
rowmax_cd34_rare <- rowMaxs(data.matrix((assays(mut_se_CD34[df$variant,])[["coverage"]] > 20) * (assays(mut_se_CD34[df$variant,])[["allele_frequency"]])))
rowmax_cd34_mgatk <- rowMaxs(data.matrix((assays(mut_se_CD34[rowData(called_mutations)$variant,])[["coverage"]] > 20) * (assays(mut_se_CD34[rowData(called_mutations)$variant,])[["allele_frequency"]])))

ggplot(data.frame(af=rowmax_cd34_rare), aes(x = af*100)) +
  geom_histogram(bins = 20, fill = "lightgrey", color = "black") +
  pretty_plot(fontsize = 7) + L_border() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("CD34")




