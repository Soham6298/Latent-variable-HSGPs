## We perform model convergence comparison for HSGPs and exact GPs for primary parameters of interest
## Need to load the output dataframe from the simulation study code before generating figures here

library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(brms)
library(parallel)
library(foreach)
library(doParallel)
#library(lemon)
#library(grid)
#library(gtable)
#library(SBC)
library(data.table)
# Import functions
source('hsgpfitfns.R')

# Import results from exact and HSGPs
compare_table <- readRDS('hsgp and exact gp sim results/hsgp_simout_se.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)
compare_x <- subset(compare_table, class == 'x' & n == 20)
compare_rho <- subset(compare_table, class == 'rho' & n == 20)
compare_alpha <- subset(compare_table, class == 'alpha' & n == 20)
compare_sigma <- subset(compare_table, class == 'sigma' & n == 20)
levels(compare_x$m)
# Change factor labels according to different simulation studies
compare_x$m <- factor(compare_x$m, levels = c('22', '26', '30', 'exact'))#c('exact', '22', '26', '30', 'pyroVI'))
compare_rho$m <- factor(compare_rho$m, levels = c('22', '26', '30', 'exact'))
compare_alpha$m <- factor(compare_alpha$m, levels = c('22', '26', '30', 'exact'))
compare_sigma$m <- factor(compare_sigma$m, levels = c('22', '26', '30', 'exact'))
# Fit summary models
formula_rhat <- bf(rhat ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

formula_bess <- bf(bess ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

formula_tess <- bf(tess ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

formula_list <- list(formula_rhat, formula_bess, formula_tess)
summary_fit <- list()
summary_fit_list <- list()
cores = 4 #50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
summary_fit = foreach(i = 1:3) %dopar% {
  library(brms)
  summary_fit_list[[i]] <- brm(formula_list[[i]], data = compare_sigma, chains = 1, cores = 1, file_refit = 'on_change')
  return(summary_fit_list)
}

m_rhat <- summary_fit[[1]][[1]]
m_bess <- summary_fit[[2]][[2]]
m_tess <- summary_fit[[3]][[3]]

#saveRDS(m_x_1000_bias, "hsgp_bias_summ_se_n1000_wideprior.rds")
#saveRDS(m_x_1000_sd, "hsgp_sd_summ_se_n1000_wideprior.rds")
## Extract summary results as conditional eff data
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

m_rhat_eff <- conditional_effects(m_rhat, effects = 'm', 
                                      conditions = make_conditions(m_rhat, 'd'),
                                      resolution = 300)
m_rhat_effs <- conditional_effects(m_rhat, effects = 'true_value:m', 
                                  conditions = make_conditions(m_rhat, 'd'),
                                  resolution = 300)
m_bess_eff <- conditional_effects(m_bess, effects = 'm', 
                                    conditions = make_conditions(m_bess, 'd'),
                                    resolution = 300)
m_bess_effs <- conditional_effects(m_bess, effects = 'true_value:m', 
                                   conditions = make_conditions(m_bess, 'd'),
                                   resolution = 300)
m_tess_eff <- conditional_effects(m_tess, effects = 'm', 
                                  conditions = make_conditions(m_tess, 'd'),
                                  resolution = 300)
m_tess_effs <- conditional_effects(m_tess, effects = 'true_value:m', 
                                   conditions = make_conditions(m_tess, 'd'),
                                   resolution = 300)

# Prepare plots
label_outdims <- c('D = 5','D = 10','D = 20')
# Check and change labels according to the number of basis functions for HSGPs
label_models <- c('HSGP(22)','HSGP(26)','HSGP(30)','Exact GP')

## Rhat
df_rhat_eff<- as.data.frame(m_rhat_eff$`m`)
levels(df_rhat_eff$cond__) <- label_outdims
levels(df_rhat_eff$effect1__) <- label_models
p_rhat_eff <- ggplot(df_rhat_eff, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Rhat') +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('(a)')

df_rhat_effs <- as.data.frame(m_rhat_effs$`true_value:m`)
levels(df_rhat_effs$cond__) <- label_outdims
levels(df_rhat_effs$effect2__) <- label_models
p_rhat_effs <- ggplot(df_rhat_effs, aes(x = effect1__, y = estimate__, 
                                          colour = effect2__, fill = effect2__)) +
  theme_bw(base_size = 35,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_rhat_effs$lower__, ymax = df_rhat_effs$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = 'Rhat', colour = 'Models', fill = 'Models') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), legend.position = 'bottom') +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))

## Bulk ESS
df_bess_eff <- as.data.frame(m_bess_eff$`m`)
levels(df_bess_eff$cond__) <- label_outdims
levels(df_bess_eff$effect1__) <- label_models
p_bess_eff <- ggplot(df_bess_eff, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Bulk-ESS') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3),
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(b)')

df_bess_effs <- as.data.frame(m_bess_effs$`true_value:m`)
levels(df_bess_effs$cond__) <- label_outdims
levels(df_bess_effs$effect2__) <- label_models
p_bess_effs <- ggplot(df_bess_effs, aes(x = effect1__, y = estimate__, 
                                        colour = effect2__, fill = effect2__)) +
  theme_bw(base_size = 35,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_bess_effs$lower__, ymax = df_bess_effs$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = 'Bulk-ESS', colour = 'Models', fill = 'Models') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), legend.position = 'bottom') +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))

## Tail ESS
df_tess_eff <- as.data.frame(m_tess_eff$`m`)
levels(df_tess_eff$cond__) <- label_outdims
levels(df_tess_eff$effect1__) <- label_models
p_tess_eff <- ggplot(df_tess_eff, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'Tail-ESS') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3),
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(c)')

df_tess_effs <- as.data.frame(m_tess_effs$`true_value:m`)
levels(df_tess_effs$cond__) <- label_outdims
levels(df_tess_effs$effect2__) <- label_models
p_tess_effs <- ggplot(df_tess_effs, aes(x = effect1__, y = estimate__, 
                                        colour = effect2__, fill = effect2__)) +
  theme_bw(base_size = 35,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_tess_effs$lower__, ymax = df_tess_effs$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = 'Tail-ESS', colour = 'Models', fill = 'Models') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), legend.position = 'bottom') +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))

# Combine the plots
p_latentx_eff <- ((p_rhat_eff + p_rhat_effs) + plot_layout(axis_titles = 'collect')) /
                    ((p_bess_eff + p_bess_effs) + plot_layout(axis_titles = 'collect'))/
                    ((p_tess_eff + p_tess_effs) + plot_layout(axis_titles = 'collect'))

ggsave('se_sigma_convcheck.pdf',
       p_latentx_eff,
       dpi = 300,
       width = 60,
       height = 60,
       units = 'cm')

