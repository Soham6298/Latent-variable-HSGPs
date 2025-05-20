## We perform model calibration for exact, HSGPs and VIGPs and check for rank based ECDF plots from simulation studies
## Need to load the output dataframe from the simulation study code before generating figures here

library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(brms)
library(lemon)
library(grid)
library(gtable)
library(SBC)
library(data.table)
# Import functions
source('hsgpfitfns.R')

# VIGP model output
gplvmpyout <- read.csv('gplvmout/GPLVMpyOut_se.csv')
summary(gplvmpyout)
colnames(gplvmpyout) <- c('n', 'd', 'sim_id','sample_id', 'true_values', 'mean', 'sd')
gplvmpyout$d <- as.factor(gplvmpyout$d)
gplvmpyout$n <- as.factor(gplvmpyout$n)
gplvmpyout$sim_id <- as.factor(gplvmpyout$sim_id)
gplvmpyout$sample_id <- as.factor(gplvmpyout$sample_id)
pyro_out <- list()
for (i in 1:nrow(gplvmpyout)){
  pyro_out[[i]] <- compare_summary_py_test(latent_mean = gplvmpyout$mean[i], 
                                           latent_sd = gplvmpyout$sd[i],
                                           true_value = gplvmpyout$true_values[i], 
                                           n_obs = gplvmpyout$n[i],  
                                           m_approx = 'pyroVI', 
                                           dims = gplvmpyout$d[i] , 
                                           variable_class = 'x', 
                                           model_name = 'pyro', 
                                           n_draws = 1000, 
                                           sim_no = gplvmpyout$sim_id[i], 
                                           sample_id = gplvmpyout$sample_id[i], 
                                           runtime = NA)
  
}
pyro_out <- rbindlist(pyro_out)
pyro_out$data_id = paste0(pyro_out$sim_id, '_', 
                          pyro_out$m, '_',
                          pyro_out$n, '_',
                          pyro_out$d)
pyro_out <- pyro_out[order(pyro_out$sim_id),]

# HSGP and exact GP model output
compare_table <- readRDS('gprout/hsgp_simout_se.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)

## Join model results for latent x
pyro_out$sim_id <- as.factor(pyro_out$sim_id)
pyro_out$n <- as.factor(pyro_out$n)
pyro_out$m <- as.factor(pyro_out$m)
pyro_out$d <- as.factor(pyro_out$d)
pyro_out$data_id <- as.factor(pyro_out$data_id)
compare_table <- rbind(compare_table, pyro_out)
compare_x <- subset(compare_table, class == 'x')
compare_x$m <- factor(compare_x$m, levels = c('exact', '22', '26', '30', 'VIGP'))

# ECDF for exact GPs
sbctestrank_exact <- compare_table[which(compare_table$n == 20 & 
                                           compare_table$d == 20 & 
                                           compare_table$class == 'x' &
                                           compare_table$m == 'exact' &
                                           compare_table$model_name == 'gp'), ]
testrank <- sbctestrank_exact[, c(1, 6, 14)]

# Rename columns
colnames(testrank) <- c('sim_id', 'variable', 'rank')

# Convert data types
testrank$rank <- as.numeric(testrank$rank)
testrank$variable <- as.factor(testrank$variable)
testrank$sim_id <- as.factor(testrank$sim_id)

# Create variable names list
x_names <- sprintf('x[%s]', seq(1:20))
plot_label <- paste0('Exact')
list_combine <- list(x_names)  
names(list_combine)[1] <- plot_label # Properly named list
# Plot for exact GP
exact_ecdf <- plot_ecdf_diff(testrank, max_rank=1000, combine_variables = list_combine) +
  theme_bw(base_size = 20, base_family = 'Times') + 
  scale_linewidth(range = 6)

# ECDF for HSGPs
sbctestrank_hsgp <- compare_table[which(compare_table$n == 20 & 
                                           compare_table$d == 20 & 
                                           compare_table$class == 'x' &
                                           compare_table$m == '22' &
                                           compare_table$model_name == 'hsgp'), ]
testrank <- sbctestrank_hsgp[, c(1, 6, 14)]

# Rename columns
colnames(testrank) <- c('sim_id', 'variable', 'rank')

# Convert data types
testrank$rank <- as.numeric(testrank$rank)
testrank$variable <- as.factor(testrank$variable)
testrank$sim_id <- as.factor(testrank$sim_id)

# Create variable names list
x_names <- sprintf('x[%s]', seq(1:20))
plot_label <- paste0('HSGP')
list_combine <- list(x_names)  
names(list_combine)[1] <- plot_label# Properly named list
# Plot for HSGPs
hsgp_ecdf <- plot_ecdf_diff(testrank, max_rank=1000, combine_variables = list_combine) +
  theme_bw(base_size = 20, base_family = 'Times') + 
  scale_linewidth(range = 6)

# ECDF for VIGPs
sbctestrank_pyro <- compare_table[which(compare_table$n == 20 & 
                                           compare_table$d == 20 & 
                                           compare_table$class == 'x' &
                                           compare_table$m == 'VIGP' &
                                           compare_table$model_name == 'pyro'), ]
testrank <- sbctestrank_pyro[, c(1, 6, 14)]

# Rename columns
colnames(testrank) <- c('sim_id', 'variable', 'rank')

# Convert data types
testrank$rank <- as.numeric(testrank$rank)
testrank$variable <- as.factor(testrank$variable)
testrank$sim_id <- as.factor(testrank$sim_id)

# Create variable names list
x_names <- sprintf('x[%s]', seq(1:20))
plot_label <- paste0('PyroVI')
list_combine <- list(x_names)  
names(list_combine)[1] <- plot_label# Properly named list
# Plot for VIGPs
pyro_ecdf <- plot_ecdf_diff(testrank, max_rank=1000, combine_variables = list_combine) +
  theme_bw(base_size = 20, base_family = 'Times') + 
  scale_linewidth(range = 6)

# Combine ECDF plots from all models
ecdf_plots <- exact_ecdf + hsgp_ecdf + pyro_ecdf + plot_layout(axes = 'collect', guides = 'collect') &
  theme(legend.position = 'bottom')

ggsave('se_ecdf_n20_d20_sbc.pdf',
       ecdf_plots,
       dpi = 300,
       width = 30,
       height = 12,
       units = 'cm')