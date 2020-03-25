library(stringr)
library(data.table)
library(dpllyr)
library(Matrix)
library(SummarizedExperiment)
load("../output/trajectory_inferences.18march2020.rda")

# Import mitomap and alter
mm <- fread("../data/mitomap_cfrm_pathogenic.tsv", header = FALSE)
mm$pos <- substring(mm[["V1"]], 2, nchar(mm[["V1"]]) - 1)
mm$var <- paste0(mm$pos, substring(mm[["V1"]], 1,1), ">", substring(mm[["V1"]], nchar(mm[["V1"]]), nchar(mm[["V1"]])))
pathogenic_variants <- as.character(mm$var)

se500 <- readRDS("../output/filtered_mitoSE_CD34-500.rds")
se800 <- readRDS("../output/filtered_mitoSE_CD34-800.rds")

# See if we overlap any pathogenic variants
intersect(rownames(se500), pathogenic_variants)
intersect(rownames(se800), pathogenic_variants)

length(intersect(rownames(se800), rownames(se500)))

# Compute raw + permuted
raw_500 <- data.frame(
  lib = str_split_fixed(rownames(colData(se500)), "-", 2)[,1],
  data.matrix(t(assays(se500)[["allele_frequency"]]))
) %>% group_by(lib) %>%
  summarize_at(vars(-group_cols()), list(mean)) %>% t()
df_500 <- data.frame(
  mutation = rownames(se500),
  pathogenic = rownames(se500) %in% pathogenic_variants,
  baseline = as.numeric(raw_500[-1,1]),
  epo  = as.numeric(raw_500[-1,2])
) %>% mutate(logFC = log2((epo + 0.000001)/(baseline+ 0.000001)))

# Permuted
set.seed(10)
perm_500 <- data.frame(
  lib = sample(str_split_fixed(rownames(colData(se500)), "-", 2)[,1]),
  data.matrix(t(assays(se500)[["allele_frequency"]]))
) %>% group_by(lib) %>%
  summarize_at(vars(-group_cols()), list(mean)) %>% t()
df_500_perm <- data.frame(
  mutation = rownames(se500),
  pathogenic = rownames(se500) %in% pathogenic_variants,
  baseline = as.numeric(perm_500[-1,1]),
  epo  = as.numeric(perm_500[-1,2])
) %>% mutate(logFC = log2((epo + 0.000001)/(baseline+ 0.000001)))


# Process the 800
raw_800 <- data.frame(
  lib = str_split_fixed(rownames(colData(se800)), "-", 2)[,1],
  data.matrix(t(assays(se800)[["allele_frequency"]]))
) %>% group_by(lib) %>%
  summarize_at(vars(-group_cols()), list(mean)) %>% t()
df_800 <- data.frame(
  mutation = rownames(se800),
  pathogenic = rownames(se800) %in% pathogenic_variants,
  baseline = as.numeric(raw_800[-1,1]),
  cc100 = as.numeric(raw_800[-1,2]),
  epo  = as.numeric(raw_800[-1,3])
) %>% mutate(logFC1 = log2((cc100 + 0.000001)/(baseline+ 0.000001)),
             logFC2 = log2((epo + 0.000001)/(baseline+ 0.000001)))


# Make trajectory plots
df_500_order <- data.frame(
  mutation = rownames(se500),
  timepoint = c(rep("A", dim(df_500)[1]), rep("B", dim(df_500)[1])),
  heteroplasmy = c(df_500$baseline, df_500$epo),
  slope = df_500$logFC >0
)

p1 <- ggplot(shuf(df_500_order), aes(x = timepoint, y = heteroplasmy, group = mutation, color = slope)) +
  geom_point(size = 0.1) + geom_line(size = 0.2) +
  scale_y_log10() + scale_x_discrete(expand = c(0.1, 0.1)) +
  scale_color_manual(values = c("firebrick", "dodgerblue4")) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "Timepoint", y = "Heteroplasmy (log scale)") +
  theme(legend.position = "none")

# Compare permuted timepoints
density_df <- data.frame(
  what = c(rep("Permuted", dim(df_500_perm)[1]), rep("Observed", dim(df_500)[1])),
  logFC = c(df_500_perm$logFC, df_500$logFC)
)

p2 <- ggplot(density_df, aes(x = logFC, fill = what)) + 
  geom_density(alpha = 0.2) +
  scale_fill_manual(values =c ("firebrick", "grey")) +
  scale_y_continuous(expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 FC", y = "empirical density") +
  theme(legend.position = "none")

ks.test(df_500_perm$logFC, df_500$logFC)


df_800$selected_var <- ifelse(df_800$mutation %in% c("3243A>T", "14322A>G", "12316G>A", "3712G>A"), df_800$mutation, "other")
p3 <- ggplot(df_800 %>% arrange(desc(selected_var)), aes(x = logFC1, y = logFC2, color = selected_var)) +
  geom_point(size = 0.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5,color = "purple") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, size = 0.5) +
  scale_y_continuous(limits = c(-5.5, 5.5)) +
  scale_x_continuous(limits = c(-5.5, 5.5)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2FC T2/T1", y = "log2FC T3/T1") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue4", "orange", "dodgerblue", "firebrick", "black"))
summary(lm(logFC2 ~ logFC1, df_800))

# Visualize trajectories of the 3 sample pool
df_800_order <- data.frame(
  mutation = rownames(se800),
  timepoint = c(rep("A", dim(df_800)[1]), rep("B", dim(df_800)[1]), rep("C", dim(df_800)[1])),
  heteroplasmy = c(df_800$baseline, df_800$cc100, df_800$epo),
  pathogenic = df_800$pathogenic, 
  selected_var = df_800$selected_var
)

# visualize individual mutation trajectories
p4 <- ggplot(shuf(df_800_order %>% filter(selected_var != "other")), aes(x = timepoint, y = heteroplasmy, group = mutation, color = mutation)) +
  geom_point(size = 0.5) + geom_line() +
  scale_y_log10() +
  labs(x = "Timepoint", y = "Heteroplasmy (log scale)") +
  pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = c("14322A>G"="dodgerblue4", "12316G>A"="firebrick", "3243A>T"="orange", "3712G>A"="dodgerblue")) +
  theme(legend.position = "none") + scale_x_discrete(expand = c(0.1, 0.2))

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, nrow = 1), 
                file = "../plots/four_plot_panel_heteroplasmy.pdf", width = 6.8, height = 1.7)


# Examine mutual variants for the supplement
l5 <- df_500 %>% filter(mutation %in%intersect(rownames(se800), rownames(se500))) %>% arrange(as.character(mutation)) %>% pull(logFC)
l8 <- df_800 %>% filter(mutation %in%intersect(rownames(se800), rownames(se500))) %>% arrange(as.character(mutation)) %>% pull(logFC1)
colonies_df <- data.frame(l5, l8)

pS <- ggplot(colonies_df, aes(x = l5, y = l8)) +
  geom_point(size = 0.2) +
  geom_smooth(method = "lm", se = FALSE, size = 0.5,color = "purple") +
  geom_abline(intercept = 0, slope = 1, linetype = 2, size = 0.5) +
  scale_y_continuous(limits = c(-4, 3)) +
  scale_x_continuous(limits = c(-4, 3)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2FC 500 input", y = "log2FC 800 input") 

summary(lm(l8 ~ l5, colonies_df))

cowplot::ggsave2(pS, file = "../plots/plot_supp_intersection_muts.pdf", width = 1.7, height = 1.7)
