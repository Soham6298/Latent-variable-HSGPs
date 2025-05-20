## We perform model evaluation for exact and HSGPs based on latent input estimation from simulation studies
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
gplvmpyout <- read.csv('gplvmout/GPLVMpyOut_se_widerho.csv')
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


# Import results from exact and HSGPs
compare_table <- readRDS('gprout/hsgp_simout_se_widerho.rds')
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
compare_table <- rbind(compare_table, pyro_out)
compare_x <- subset(compare_table, class == 'x')
compare_rho <- subset(compare_table, class == 'rho')
compare_alpha <- subset(compare_table, class == 'alpha')
compare_sigma <- subset(compare_table, class == 'sigma')
levels(compare_x$m)
compare_x$m <- factor(compare_x$m, levels = c('exact', '22', '26', '30', 'pyroVI'))
compare_x_20 <- subset(compare_x, n == 20)
compare_x_50 <- subset(compare_x, n == 50)
compare_x_200 <- subset(compare_x, n == 200)

# Fit summary models
formula_bias <- bf(abs_bias ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
                   sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

formula_sd <- bf(sd ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
                 sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

m_x_20_bias <- brm(formula_bias, data = compare_x_20, chains = 2, cores = 2, file_refit = 'on_change')

m_x_20_sd <- brm(formula_sd,data = compare_x_20, chains = 2, cores = 2, file_refit = 'on_change')

m_x_50_bias <- brm(formula_bias,data = compare_x_50, chains = 2, cores = 2, file_refit = 'on_change')

m_x_50_sd <- brm(formula_sd, data = compare_x_50, chains = 2, cores = 2, file_refit = 'on_change')

m_x_200_bias <- brm(formula_bias,data = compare_x_200, chains = 2, cores = 2, file_refit = 'on_change')

m_x_200_sd <- brm(formula_sd,data = compare_x_200, chains = 2, cores = 2, file_refit = 'on_change')


## Extract summary results as conditional eff data
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

m_bias_eff20 <- conditional_effects(m_x_20_bias, effects = 'm', 
                                    conditions = make_conditions(m_x_20_bias, 'd'),
                                    resolution = 300)
m_bias_eff20_s <- conditional_effects(m_x_20_bias, effects = 'true_value:m', 
                                      conditions = make_conditions(m_x_20_bias, 'd'),
                                      resolution = 300)
m_sd_eff20 <- conditional_effects(m_x_20_sd, effects = 'm', 
                                  conditions = make_conditions(m_x_20_sd, 'd'),
                                  resolution = 300)
m_sd_eff20_s <- conditional_effects(m_x_20_sd, effects = 'true_value:m', 
                                    conditions = make_conditions(m_x_20_sd, 'd'),
                                    resolution = 300)
m_bias_eff50 <- conditional_effects(m_x_50_bias, effects = 'm', 
                                    conditions = make_conditions(m_x_50_bias, 'd'),
                                    resolution = 300)
m_bias_eff50_s <- conditional_effects(m_x_50_bias, effects = 'true_value:m', 
                                      conditions = make_conditions(m_x_50_bias, 'd'),
                                      resolution = 300)
m_sd_eff50 <- conditional_effects(m_x_50_sd, effects = 'm', 
                                  conditions = make_conditions(m_x_50_sd, 'd'),
                                  resolution = 300)
m_sd_eff50_s <- conditional_effects(m_x_50_sd, effects = 'true_value:m', 
                                    conditions = make_conditions(m_x_50_sd, 'd'),
                                    resolution = 300)
m_bias_eff200 <- conditional_effects(m_x_200_bias, effects = 'm', 
                                     conditions = make_conditions(m_x_200_bias, 'd'),
                                     resolution = 300)
m_bias_eff200_s <- conditional_effects(m_x_200_bias, effects = 'true_value:m', 
                                       conditions = make_conditions(m_x_200_bias, 'd'),
                                       resolution = 300)
m_sd_eff200 <- conditional_effects(m_x_200_sd, effects = 'm', 
                                   conditions = make_conditions(m_x_200_sd, 'd'),
                                   resolution = 300)
m_sd_eff200_s <- conditional_effects(m_x_200_sd, effects = 'true_value:m', 
                                     conditions = make_conditions(m_x_200_sd, 'd'),
                                     resolution = 300)

# Prepare plots
label_outdims <- c('D = 5','D = 10','D = 20')
label_models <- c('Exact', 'HSGP(22)','HSGP(26)','HSGP(30)','VIGP')

# Posterior bias plots
df_bias_eff20 <- as.data.frame(m_bias_eff20$`m`)
levels(df_bias_eff20$cond__) <- label_outdims
levels(df_bias_eff20$effect1__) <- label_models
p_bias_eff20 <- ggplot(df_bias_eff20, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=25,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Abs Bias') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('(a) N = 20')

df_bias_eff50 <- as.data.frame(m_bias_eff50$`m`)
levels(df_bias_eff50$cond__) <- label_outdims
levels(df_bias_eff50$effect1__) <- label_models[-1]
p_bias_eff50 <- ggplot(df_bias_eff50, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=25,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Abs Bias') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('(b) N = 50')

df_bias_eff200 <- as.data.frame(m_bias_eff200$`m`)
levels(df_bias_eff200$cond__) <- label_outdims
levels(df_bias_eff200$effect1__) <- label_models[-1]
p_bias_eff200 <- ggplot(df_bias_eff200, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=25,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Abs Bias') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('(c) N = 200')

# Posterior SD plots
df_sd_eff20 <- as.data.frame(m_sd_eff20$`m`)
levels(df_sd_eff20$cond__) <- label_outdims
levels(df_sd_eff20$effect1__) <- label_models
p_sd_eff20 <- ggplot(df_sd_eff20, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=25,
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
  scale_colour_manual(values = c("#000000" ))

df_sd_eff50 <- as.data.frame(m_sd_eff50$`m`)
levels(df_sd_eff50$cond__) <- label_outdims
levels(df_sd_eff50$effect1__) <- label_models[-1]
p_sd_eff50 <- ggplot(df_sd_eff50, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=25,
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
  scale_colour_manual(values = c("#000000" ))

df_sd_eff200 <- as.data.frame(m_sd_eff200$`m`)
levels(df_sd_eff200$cond__) <- label_outdims
levels(df_sd_eff200$effect1__) <- label_models[-1]
p_sd_eff200 <- ggplot(df_sd_eff200, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=25,
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
p_latentx_eff <- ((p_bias_eff20 + p_sd_eff20) / 
                      (p_bias_eff50 + p_sd_eff50) / 
                      (p_bias_eff200 + p_sd_eff200)) + 
                       plot_layout(axis_titles = 'collect')

ggsave('se_widerho_latentx.pdf',
       p_latentx_eff,
       dpi = 300,
       width = 60,
       height = 40,
       units = 'cm')
