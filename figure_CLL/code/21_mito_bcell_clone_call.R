library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(BuenColors)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(ggrepel)
library(Seurat)

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)
pt1 <- "../CLL_CD19pos_CR-mtMask_mgatk/final/CLL_CD19pos_CR-mtMask_mgatk.rds"
pt2 <- "../Patient2_CLL_CD19_pos_CR-mtMask_mgatk/final/Patient2_CLL_CD19_pos_CR-mtMask_mgatk.rds"
SE <- readRDS(pt2)
SE <- SE[,colData(SE)$depth >= 20]
dim(SE)

call_mutations_mgatk <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10){
  
  # Determinie key coverage statistics every which way
  cov <- assays(SE)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  
  # Process mutation for one alternate letter
  process_letter <- function(letter){
    print(letter)
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]
    
    # Single cell functions
    getMutMatrix <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_fw  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]]) / cov_fw)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_rev  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_rev")]]) / cov_rev)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    # Bulk functions
    getBulk <- function(letter){
      vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / Matrix::rowSums(cov))[boo]
      return(vec)
    }
    rowVars <- function(x, ...) {
      Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    }
    
    update_missing_w_zero <- function(vec){
      ifelse(is.na(vec)  | is.nan(vec), 0, vec)
    }
    # Set up correlation per non-zero mutation based on the strands
    dt <- merge(data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_fw")]][boo,])), 
                data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_rev")]][boo,])), 
                by.x = c("i", "j"), by.y = c("i", "j"), 
                all = TRUE)[x.x >0 | x.y >0]
    dt$x.x <- update_missing_w_zero(dt$x.x)
    dt$x.y <- update_missing_w_zero(dt$x.y)
    
    dt2 <- data.table(variant = variant_name[dt[[1]]],
                      cell_idx = dt[[2]], 
                      forward = dt[[3]],
                      reverse = dt[[4]])
    rm(dt)
    cor_dt <- dt2[, .(cor = cor(c(forward), c(reverse), method = "pearson", use = "pairwise.complete")), by = list(variant)]
    
    # Put in vector for convenience
    cor_vec_val <- cor_dt$cor
    names(cor_vec_val) <- as.character(cor_dt$variant )
    
    # Compute the single-cell data
    mat <- getMutMatrix(letter)
    mmat <- sparseMatrix(
      i = c(summary(mat)$i,dim(mat)[1]),
      j = c(summary(mat)$j,dim(mat)[2]),
      x = c(update_missing_w_zero(summary(mat)$x), 0)
    )
    
    # Compute bulk statistics
    mean = update_missing_w_zero(getBulk(letter))
    
    # Stablize variances by replacing low coverage cells with mean
    if(stabilize_variance){
      
      # Get indices of cell/variants where the coverage is low and pull the mean for that variant
      idx_mat <- which(data.matrix(cov[boo,] < low_coverage_threshold), arr.ind = TRUE)
      idx_mat_mean <- mean[idx_mat[,1]]
      
      # Now, make sparse matrices for quick conversion
      ones <- 1 - sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = 1
      )
      
      means_mat <- sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = c(idx_mat_mean, 0)
      )
      
      mmat2 <- mmat * ones + means_mat
      variance = rowVars(mmat2)
      rm(mmat2); rm(ones); rm(means_mat); rm(idx_mat); rm(idx_mat_mean)
      
    } else {
      variance = rowVars(mmat)
    }
    
    detected <- (assays(SE)[[paste0(letter, "_counts_fw")]][boo,] >= 2) + (assays(SE)[[paste0(letter, "_counts_rev")]][boo,] >=2 )
    
    # Compute per-mutation summary statistics
    var_summary_df <- data.frame(
      position = position_filt,
      nucleotide = nucleotide, 
      variant = variant_name,
      vmr = variance/(mean + 0.00000000001),
      mean = round(mean,7),
      variance = round(variance,7),
      n_cells_detected = Matrix::rowSums(detected == 2),
      n_cells_over_5 = Matrix::rowSums(mmat >= 0.05), 
      n_cells_over_10 = Matrix::rowSums(mmat >= 0.10),
      n_cells_over_20 = Matrix::rowSums(mmat >= 0.20),
      strand_correlation = cor_vec_val[variant_name],
      mean_coverage = Matrix::rowMeans(cov)[boo], 
      stringsAsFactors = FALSE, row.names = variant_name
    )
    se_new <- SummarizedExperiment(
      rowData = var_summary_df, 
      colData = colData(SE), 
      assays = list(allele_frequency = mmat, coverage = cov[boo,])
    )
    return(se_new)
  }
  
  return(SummarizedExperiment::rbind(process_letter("A"), process_letter("C"), process_letter("G"), process_letter("T")))
  
}

mut_se <- call_mutations_mgatk(SE, stabilize_variance = FALSE)
misc_df <- data.frame(rowData(mut_se))


filter_df <- misc_df %>%  filter(n_cells_detected >= 5 & strand_correlation >= 0.65 & log10(vmr) > -2)
dim(filter_df)

p1 <- ggplot(misc_df %>% filter(n_cells_detected >= 5), aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.6)) +
  geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
  labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
  pretty_plot(fontsize = 7) + L_border() +
  theme(legend.position = "bottom") + 
  geom_vline(xintercept = 0.6, linetype = 2) +
  geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
cowplot::ggsave(p1, file = "output/CLL_Pt2_CD19pos.pdf", width = 1.2, height = 1.2)

afp <- as.matrix(data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),]))
afp[afp > 0.1] <- 0.1

# Infer clusters
seuratSNN <- function(matSVD, resolution = 0.8,k.param = 20){
  set.seed(1)
  rownames(matSVD) <- make.unique(rownames(matSVD))
  obj <- FindNeighbors(matSVD, k.param = k.param)
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

clusters <- seuratSNN(t(sqrt(afp)), resolution = 1.2, k.param = 50) # default for PT1; resolution = 1.2; k.param = 50 for PT2
table(clusters)
clusters <- as.character(as.numeric(factor(clusters)))
re_map <- rev(c("a", "b", "c", "d", "e", "f")); names(re_map) <- c("1", "2", "3", "5", "6","4")
clusters <- re_map[clusters]

# Determine allele frequencies per cell
average_pretty <- function(vec) return(round(sum(vec)/length(vec) * 100, 3))
mut_df <- data.frame(
  clusters, 
  t(as.matrix(data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),])))
) %>% group_by(clusters)%>%
  summarise_at(vars(-group_cols("clusters")), average_pretty) %>% data.frame

# Create pretty colors
set.seed(2)
names_clusters <- unique(sort((clusters)))
vec_go <- c("grey", jdb_palettes[["flame_light"]], "dodgerblue2")[1:length(unique(clusters))]
names(vec_go) <- names_clusters

# Make data.frame for stuff
df <- data.frame(
  cell_id = colnames(afp), 
  cluster_id = clusters
) %>% arrange(desc(clusters))

ha_col <- HeatmapAnnotation(cell = as.character(df$cluster_id),
                            col = list(cell = vec_go))
#afp[afp < 0.01] <- 0 

var_order_1 <- c("12213G>A", "4716C>A", "2587G>A", "5274G>A", "3179G>A", "2702G>A", "11711G>A", 
                 "9950T>C","5140G>A",  "16519T>C", "16052C>T","5809G>A",
                 "14858G>A", "10191T>C", "5835G>A", "1872T>C", "16390G>A", "1260A>G", "12067C>T"
                 )

var_order_2 <- c("4853G>A","12980G>A",  "7775G>A", 
                 "3496G>A", "4429G>A", "4972G>A",
                 "2702G>A", "14641A>G",
                 "1031G>A","1782G>A","2989G>A", "3244G>A" , "5293G>A" , "8020G>A" ,"8858G>A" , "11711G>A" ,"12236G>A" ,"13513G>A" ,"15437G>A" ,"16390G>A", "16391G>A" ,"2275T>C" , "11982T>C" ,"10367C>T")
setdiff(var_order_2, rownames(afp))
setdiff(rownames(afp), var_order_2)
var_order <- var_order_2

# Make heatmap
png(paste0("output/PT2.png"), width=4, height=2, unit = "in", res = 500)
hm <- Heatmap((data.matrix(afp[var_order,as.character(df$cell_id)])), #
              col=as.character(jdb_palette("solar_rojos",type="continuous")),
              show_row_names = TRUE, 
              top_annotation=ha_col,
              cluster_columns = FALSE,
              name = "AF",
              row_names_gp = gpar(fontsize = 4),
              cluster_rows = FALSE, 
              show_column_names = FALSE)
hm 
dev.off()

pdf(paste0("output/PT2.pdf"), width=4, height=2)
hm 
dev.off()

saveRDS(df, file = "output/PT2_clone_definition.rds")


if(FALSE){
  mut_df_o <- data.frame(
    cluster = clusters, 
    t(data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),]))
  ) 
  saveRDS(mut_df_o, file = "output/PT1_specificVariants_forSupplement.rds")
  ggplot(mut_df_o, aes(x = cluster, y = X12067C.T)) + geom_quasirandom()
  
}

if(FALSE){
  mut_df_o <- data.frame(
    cluster = clusters, 
    t(data.matrix(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant),]))
  ) 
  saveRDS(mut_df_o, file = "output/PT2_specificVariants_forSupplement.rds")
  ggplot(mut_df_o, aes(x = cluster, y = X10367C.T)) + geom_quasirandom()
  
}
