## We perform model evaluation for HSGPs based on latent input estimation from simulation studies with N = 1000
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

# Import results from VIGPs
gplvmpyout <- read.csv('vigp sim results/GPLVMpyOut_se_n1000.csv')
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

gplvmpyout_fixpar <- read.csv('vigp sim results/GPLVMpyOut_se_fixedparams_n1000.csv')
colnames(gplvmpyout_fixpar) <- c('n', 'd', 'sim_id','sample_id', 'true_values', 'mean', 'sd')
gplvmpyout_fixpar$d <- as.factor(gplvmpyout_fixpar$d)
gplvmpyout_fixpar$n <- as.factor(gplvmpyout_fixpar$n)
gplvmpyout_fixpar$sim_id <- as.factor(gplvmpyout_fixpar$sim_id)
gplvmpyout_fixpar$sample_id <- as.factor(gplvmpyout_fixpar$sample_id)
pyro_out_fixpar <- list()
for (i in 1:nrow(gplvmpyout)){
  pyro_out_fixpar[[i]] <- compare_summary_py_test(latent_mean = gplvmpyout_fixpar$mean[i], 
                                           latent_sd = gplvmpyout_fixpar$sd[i],
                                           true_value = gplvmpyout_fixpar$true_values[i], 
                                           n_obs = gplvmpyout_fixpar$n[i],  
                                           m_approx = 'pyroVI_fix', 
                                           dims = gplvmpyout_fixpar$d[i] , 
                                           variable_class = 'x', 
                                           model_name = 'pyro', 
                                           n_draws = 1000, 
                                           sim_no = gplvmpyout_fixpar$sim_id[i], 
                                           sample_id = gplvmpyout_fixpar$sample_id[i], 
                                           runtime = NA)
  
}
pyro_out_fixpar <- rbindlist(pyro_out_fixpar)
pyro_out_fixpar$data_id = paste0(pyro_out_fixpar$sim_id, '_', 
                          pyro_out_fixpar$m, '_',
                          pyro_out_fixpar$n, '_',
                          pyro_out_fixpar$d)
pyro_out_fixpar <- pyro_out_fixpar[order(pyro_out_fixpar$sim_id),]


# Import results from exact and HSGPs
compare_table <- readRDS('hsgp and exact gp sim results/hsgp_simout_se_n1000_wideprior.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)
## Join all model results for latent x
pyro_out$sim_id <- as.factor(pyro_out$sim_id)
pyro_out$n <- as.factor(pyro_out$n)
pyro_out$m <- as.factor(pyro_out$m)
pyro_out$d <- as.factor(pyro_out$d)
pyro_out$data_id <- as.factor(pyro_out$data_id)
pyro_out_fixpar$sim_id <- as.factor(pyro_out_fixpar$sim_id)
pyro_out_fixpar$n <- as.factor(pyro_out_fixpar$n)
pyro_out_fixpar$m <- as.factor(pyro_out_fixpar$m)
pyro_out_fixpar$d <- as.factor(pyro_out_fixpar$d)
pyro_out_fixpar$data_id <- as.factor(pyro_out_fixpar$data_id)
compare_table <- rbind(compare_table, pyro_out, pyro_out_fixpar)
compare_x <- subset(compare_table, class == 'x')
compare_rho <- subset(compare_table, class == 'rho')
compare_alpha <- subset(compare_table, class == 'alpha')
compare_sigma <- subset(compare_table, class == 'sigma')
levels(compare_x$m)
# Change factor labels according to different simulation studies
compare_x$m <- factor(compare_x$m, levels = c('15', '19', '23', 'pyroVI', 'pyroVI_fix'))#c('exact', '22', '26', '30', 'pyroVI'))

# Fit summary models
formula_bias <- bf(abs_bias ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
                   sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

formula_sd <- bf(sd ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
                 sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))


m_x_1000_bias <- brm(formula_bias,data = compare_x, chains = 2, cores = 2, file_refit = 'on_change')
m_x_1000_sd <- brm(formula_sd,data = compare_x, chains = 2, cores = 2, file_refit = 'on_change')

#saveRDS(m_x_1000_bias, "hsgp_bias_summ_se_n1000_wideprior.rds")
#saveRDS(m_x_1000_sd, "hsgp_sd_summ_se_n1000_wideprior.rds")
## Extract summary results as conditional eff data
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

m_bias_eff1000 <- conditional_effects(m_x_1000_bias, effects = 'm', 
                                      conditions = make_conditions(m_x_1000_bias, 'd'),
                                      resolution = 300)
m_sd_eff1000 <- conditional_effects(m_x_1000_sd, effects = 'm', 
                                    conditions = make_conditions(m_x_1000_sd, 'd'),
                                    resolution = 300)

# Prepare plots
label_outdims <- c('D = 5','D = 10','D = 20')
# Check and change labels according to the number of basis functions for HSGPs
label_models <- c('HSGP(15)','HSGP(19)','HSGP(23)','VIGP')

df_bias_eff1000 <- as.data.frame(m_bias_eff1000$`m`)
levels(df_bias_eff1000$cond__) <- label_outdims
levels(df_bias_eff1000$effect1__) <- label_models
p_bias_eff1000 <- ggplot(df_bias_eff1000, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Bias') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('N = 1000')

df_sd_eff1000 <- as.data.frame(m_sd_eff1000$`m`)
levels(df_sd_eff1000$cond__) <- label_outdims
levels(df_sd_eff1000$effect1__) <- label_models
p_sd_eff1000 <- ggplot(df_sd_eff1000, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'SD') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3),
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000"))

# Combine the plots
p_latentx_eff <- (p_bias_eff1000 + p_sd_eff1000) + plot_layout(axis_titles = 'collect')

ggsave('se_latentx_n1000_wideprior.pdf',
       p_latentx_eff,
       dpi = 300,
       width = 60,
       height = 20,
       units = 'cm')

