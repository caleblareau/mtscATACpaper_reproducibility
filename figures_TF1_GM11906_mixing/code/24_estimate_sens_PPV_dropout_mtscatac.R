library(dplyr)
library(data.table)
library(BuenColors)

# Global parameters
contamination_percent <- 0.19 # essentially our noise parameter estimated from the data
n_cells_population <- 1000 # total cells
n_cells_subclone <- 100 # subset of cells with the mutation
n_cells_background <- n_cells_population - n_cells_subclone

# Function to estimate Sensitivity and PPV for one mutation
estimate_sensitity_and_ppv <- function(heteroplasmy, coverage, n_vals = 10000){
  
  # Compute detection threshold
  threshold = (heteroplasmy / 2 / 100)
  
  set.seed(123 + coverage)
  
  # Simulate data based on putative contamination rates
  
  # Foreground is number of times we observe the mutation per cell in clone... binarize for detection
  foreground <- sapply(1:n_vals, function(x) rbinom(n_cells_subclone, coverage, (heteroplasmy)/100))/coverage > threshold
  
  # Background is number of times we errantly observe the mutation in the cells that weren't specificied to have the mutation
  background <- sapply(1:n_vals, function(x) rbinom(n_cells_background, coverage, (contamination_percent)/100))/coverage > threshold
  
  # Compute per simulation
  sens_vec <- colMeans(foreground)
  ppv_vec <- colSums(foreground)/(colSums(foreground) + colSums(background))
  
  # Export values
  data.frame(
    sensitivity = round(mean(sens_vec),4),
    PPV = round(mean(ppv_vec), 4),
    sensitivity_down = round(quantile(sens_vec, 0.025),4),
    sensitivity_up = round(quantile(sens_vec, 0.975),4),
    PPV_down = round(quantile(sens_vec, 0.025), 4),
    PPV_up = round(quantile(ppv_vec, 0.975),4),
    heteroplasmy,
    coverage
  )
}


# Second iteration-- at least 5% heteroplasmy
lapply(c(20, 50, 100), function(coverage){
  lapply(c(2, 5, 15, 25, 35, 45), function(heteroplasmy){
    estimate_sensitity_and_ppv(heteroplasmy,coverage)
  }) %>% rbindlist() %>% data.frame()
})%>% rbindlist() %>% data.frame() -> all_df

p1 <- ggplot(all_df, aes(x = heteroplasmy, y = sensitivity, color = as.factor(coverage))) +
  geom_line() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + pretty_plot(fontsize = 8) +
  L_border() + labs(x = "% Heteroplasmy", y = "Sensitivity", color = "cov") +
  theme(legend.position = "bottom") +
  scale_color_manual(values = jdb_palette("brewer_spectra")[c(1,2,9)]) +
  geom_vline(xintercept = c(5, 15), linetype = 2) 

p2 <- ggplot(all_df, aes(x = heteroplasmy, y = PPV, color = as.factor(coverage))) +
  geom_line() + scale_y_continuous(limits = c(0,1), expand = c(0,0)) + pretty_plot(fontsize = 8) +
  L_border() + labs(x = "% Heteroplasmy", y = "Positive Predictive Value", color = "cov") +
  theme(legend.position = "bottom") + scale_color_manual(values = jdb_palette("brewer_spectra")[c(1,2,9)]) +
  geom_vline(xintercept = c(5, 15), linetype = 2)

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1), file = "../plots/PPV_Sens_simulation.pdf", width = 3.6, height = 2.2)

all_df %>% filter(heteroplasmy %in% c(5,15))


#-------
# Now look at dropout
#-------

estimate_dropout <- function(heteroplasmy, coverage, n_vals = 10000){
  set.seed(123 + coverage)
  
  # Dropout is the number of times that we fully lose the allele based on coverage
  dropout <- mean(rbinom(n_vals, coverage, (heteroplasmy)/100) == 0)
  data.frame(
    dropout,
    heteroplasmy,
    coverage
  )
}


# Loop over the same values for dropout
lapply(c(20, 50, 100), function(coverage){
  lapply(c(1:45), function(heteroplasmy){
    estimate_dropout(heteroplasmy,coverage)
  }) %>% rbindlist() %>% data.frame()
})%>% rbindlist() %>% data.frame() -> d_df


cpD <- ggplot(d_df, aes(x = heteroplasmy, y = dropout*100, color = as.factor(coverage))) +
  geom_line() + scale_y_continuous(limits = c(0,100), expand = c(0,0)) + pretty_plot(fontsize = 8) +
  L_border() + labs(x = "% Heteroplasmy", y = "% dropout", color = "cov") +
  theme(legend.position = "bottom") + scale_color_manual(values = jdb_palette("brewer_spectra")[c(1,2,9)]) +
  geom_vline(xintercept = c(5), linetype = 2)

cowplot::ggsave2(cpD, file = "../plots/dropout.pdf", width = 1.8, height = 2.2)

