library(dplyr)
library(data.table)
library(stringr)

# Variants for each
tf1_vars <- c("9824T>C", "6680T>C", "8701A>G", "15301G>A", "12811T>C", "14783T>C", "4164A>G", "7853G>A", 
              "7684T>C", "10345T>C", "12705C>T", "15043G>A", "5351A>G", "10873T>C", "150C>T", "6455C>T",
              "9540T>C", "12405C>T", "199T>C", "16129G>A", "10400C>T", "8552T>C", "4071C>T", "4048G>A",
              "16298T>C", "16297T>C", "16223C>T", "5460G>A", "16189T>C")
gm11906_vars <- c("6734G>A", "3010G>A", "6554C>T", "14798T>C", "11251A>G", "15452C>A", "12612A>G", "16069C>T",
                  "16092T>C", "16126T>C", "12127G>A", "462C>T", "4216T>C", "16261C>T", "228G>A", "295C>T")

# Function to parse out the parts of the string as necessary
substrRight <- function(x, n =1){
  substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n = 1){
  substr(x, 1, nchar(x)-n)
}

# Pull the letter for each cell line
# For TF1 variants, its pos_MERRF>TF1
tf1temp <- str_split_fixed(tf1_vars, ">", 2)
tf1var_df <- data.frame(
  pos = as.numeric(substrLeft(tf1temp[,1])),
  MERRF = substrRight(tf1temp[,1]),
  TF1 = tf1temp[,2]
)
# For MERRF variants, its pos_TF1>MERRF
mtemp <- str_split_fixed(gm11906_vars, ">", 2)
mvar_df <- data.frame(
  pos = as.numeric(substrLeft(mtemp[,1])),
  TF1 = substrRight(mtemp[,1]),
  MERRF = mtemp[,2]
)

# Define a global variable of all of the variants 
df2 <- rbind(mvar_df, tf1var_df)

# Function to get the per barcode counts per cell type (MERRF/TF1) by alt letter
extractme <- function(cell, letter, SE){
  idxx <- df2[which(df2[,cell] == letter), "pos"]
  Matrix::colSums(assays(SE)[[paste0(letter, "_counts_fw")]][idxx, ] + assays(SE)[[paste0(letter, "_counts_rev")]][idxx, ])
}

# Given a summarized experiment from mgatk, compute the essentials for ultimately determining contamination
process_SE_contamination <- function(SE, library){
  df3 <- data.frame(
    barcode = colnames(SE),
    TF1 = extractme("TF1", "A", SE) +  extractme("TF1", "C", SE) +  extractme("TF1", "G", SE) +  extractme("TF1", "T", SE),
    MERRF = extractme("MERRF", "A", SE) +  extractme("MERRF", "C", SE) +  extractme("MERRF", "G", SE) +  extractme("MERRF", "T", SE), 
    depth = colData(SE)$depth
  )
  df3 <-  df3 %>% mutate(minor_population = pmin(TF1/(MERRF + TF1 + 0.001)*100 ,MERRF/(MERRF + TF1 + 0.001)*100), library = library)
  df3                        
}

# Function to estimate contamination based on reads
estimate_contamination <- function(df3, rm_doublets = FALSE){
  if(rm_doublets){
    df3 <- df3[df3$minor_population < 5,]
  }
  sum(pmin(df3$TF1, df3$MERRF)) / sum(df3$TF1 + df3$MERRF)*100
}

noFA <- process_SE_contamination(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_2_NP40_CR-mtMask_mgatk.rds"), "noFA")
noFA$density <- get_density(log10(noFA$TF1 + 1), log10(noFA$MERRF + 1))

ggplot(noFA %>% filter(TF1 > 0 & MERRF > 0), aes(x = TF1, y = MERRF, color = minor_population)) +
  geom_point() + scale_y_log10(breaks = c(10, 100, 1000)) + scale_x_log10(breaks = c(10, 100, 1000)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_gradientn(colors = jdb_palette("flame_light")) +
  theme(legend.position = "none")

FA <- rbind(
  process_SE_contamination(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_1h_CR-mtMask_mgatk.rds"), "FA"),
  process_SE_contamination(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_Fix_6h_CR-mtMask_mgatk.rds"), "FA")
)

estimate_contamination(noFA, FALSE)
estimate_contamination(FA, TRUE)
estimate_contamination(FA, FALSE)


ggplot(FA %>% filter(TF1 > 0 & MERRF > 0), aes(x = TF1, y = MERRF, color = minor_population)) +
  geom_point() + scale_y_log10(breaks = c(10, 100, 1000), expand = c(0,0)) + scale_x_log10(breaks = c(10, 100, 1000), expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  scale_color_gradientn(colors = jdb_palette("flame_light")) +
  theme(legend.position = "none")


#----
# Examine the variable fixation conditions
#----

fix1 <- process_SE_contamination(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_3_Fix_12h_CR-mtMask_mgatk.rds"), "fix1")
fix01 <- process_SE_contamination(readRDS("../../../mtscATACpaper_large_data_files/source/mgatk_output/Mix_4_Fix_12h_CR-mtMask_mgatk.rds"), "fix01")

estimate_contamination(fix1, FALSE)
estimate_contamination(fix01, FALSE)
