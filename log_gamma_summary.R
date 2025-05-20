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

# Import results from exact and HSGPs
compare_table <- readRDS('gprout/hsgp_simout_se.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)

# Combine results for latent x
pyro_out$sim_id <- as.factor(pyro_out$sim_id)
pyro_out$n <- as.factor(pyro_out$n)
pyro_out$m <- as.factor(pyro_out$m)
pyro_out$d <- as.factor(pyro_out$d)
pyro_out$data_id <- as.factor(pyro_out$data_id)
compare_table <- rbind(compare_table, pyro_out)
compare_x <- subset(compare_table, class == 'x')
compare_x$m <- factor(compare_x$m, levels = c('exact', '22', '26', '30', 'pyroVI'))
# Create indicator for each simulation conditions
compare_x$subset_id <- paste0(compare_x$n,'_',compare_x$m,'_',compare_x$d)
# Separate results accordingly
list_df <- split(compare_x, compare_x$subset_id)
# Compute log gamma scores for each condition
gamma_stat <- list()
for(i in 1:39){
  if(list_df[[i]]$n[1] == '20'){
    gamma_stat[[i]] <- rep(NA, 20)
    for(j in 1:20){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('x','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
  }else if(list_df[[i]]$n[1] == '50'){
    gamma_stat[[i]] <- rep(NA, 50)
    for(j in 1:50){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('x','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
  }else if(list_df[[i]]$n[1] == '200'){
    gamma_stat[[i]] <- rep(NA, 200)
    for(j in 1:200){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('x','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
  }
}
# Create long table formate for the log gamma scores
gamma_stat_long <- list()
for(i in 1:39){
  if(list_df[[i]]$n[1] == '20'){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
  }else if(list_df[[i]]$n[1] == '50'){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
  }else if(list_df[[i]]$n[1] == '200'){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
  }
}
# Combine the log gamma scores
log_gamma_stat <- rbindlist(gamma_stat_long)
# Clean the results (if required)
log_gamma_stat$log_gamma = ifelse(is.infinite(log_gamma_stat$log_gamma),
                                  -1000, log_gamma_stat$log_gamma)
# Compute log gamma - critical value
log_gamma_stat$log_gamma_diff <- log_gamma_stat$log_gamma - log(SBC:::adjust_gamma(N = 50, K = 1000, L = 1))
# Subset the results for unique values
log_gamma_data <- subset(log_gamma_stat, sim_id==1) 
# Separate log gamma scores per sample size
log_gamma_data_20 <- subset(log_gamma_data, n == 20)
log_gamma_data_50 <- subset(log_gamma_data, n == 50)
log_gamma_data_200 <- subset(log_gamma_data, n == 200)

# Fit summary models to log gamma scores
m_log_gamma20 <- brm(bf(log_gamma_diff ~ (1 + m) * d + s(true_value, by = m),
                        sigma ~ (1 + m) * d + s(true_value, by = m)),
                     data = log_gamma_data_20, chains = 2, cores = 2, file_refit = 'on_change')

m_log_gamma50 <- brm(bf(log_gamma_diff ~ (1 + m) * d + s(true_value, by = m),
                        sigma ~ (1 + m) * d + s(true_value, by = m)),
                     data = log_gamma_data_50, chains = 2, cores = 2, file_refit = 'on_change')

m_log_gamma200 <- brm(bf(log_gamma_diff ~ (1 + m) * d + s(true_value, by = m),
                         sigma ~ (1 + m) * d + s(true_value, by = m)),
                      data = log_gamma_data_200, chains = 2, cores = 2, file_refit = 'on_change')

# Extract results as conditional effects
m_log_gamma_eff20 <- conditional_effects(m_log_gamma20, effects = 'm', 
                                         conditions = make_conditions(m_log_gamma20, 'd'),
                                         resolution = 300)
m_log_gamma_eff50 <- conditional_effects(m_log_gamma50, effects = 'm', 
                                         conditions = make_conditions(m_log_gamma50, 'd'),
                                         resolution = 300)
m_log_gamma_eff200 <- conditional_effects(m_log_gamma200, effects = 'm', 
                                          conditions = make_conditions(m_log_gamma200, 'd'),
                                          resolution = 300)
m_log_gamma_eff20_s <- conditional_effects(m_log_gamma20, effects = 'true_value:m', 
                                           conditions = make_conditions(m_log_gamma20, 'd'),
                                           resolution = 300)#, method = 'posterior_predict')
m_log_gamma_eff50_s <- conditional_effects(m_log_gamma50, effects = 'true_value:m', 
                                           conditions = make_conditions(m_log_gamma50, 'd'),
                                           resolution = 300)#,  method = 'posterior_predict')
m_log_gamma_eff200_s <- conditional_effects(m_log_gamma200, effects = 'true_value:m', 
                                            conditions = make_conditions(m_log_gamma200, 'd'),
                                            resolution = 300)#,  method = 'posterior_predict')

cols <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# Generate plots
# errorbar plots
label_outdims <- c('D = 5','D = 10','D = 20')
label_models <- c('Exact', 'HSGP(22)','HSGP(26)','HSGP(30)','VIGP')
df_log_gamma_eff20 <- as.data.frame(m_log_gamma_eff20$m)
levels(df_log_gamma_eff20$cond__) <- label_outdims
levels(df_log_gamma_eff20$effect1__) <- label_models
# Plot for N = 20
p_log_gamma_eff20 <- ggplot(df_log_gamma_eff20, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~(differenced))) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(a) N = 20')

# Plot for N = 50
df_log_gamma_eff50 <- as.data.frame(m_log_gamma_eff50$m)
levels(df_log_gamma_eff50$cond__) <- label_outdims
levels(df_log_gamma_eff50$effect1__) <- label_models[-1]
p_log_gamma_eff50 <- ggplot(df_log_gamma_eff50, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~(differenced))) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6), axis.title.x = element_blank()) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(b) N = 50')

# Plot for N = 200
df_log_gamma_eff200 <- as.data.frame(m_log_gamma_eff200$m)
levels(df_log_gamma_eff200$cond__) <- label_outdims
levels(df_log_gamma_eff200$effect1__) <- label_models[-1]
p_log_gamma_eff200 <- ggplot(df_log_gamma_eff200, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~(differenced))) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6),
        axis.title.y = element_blank()) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(c) N = 200')

# Spline plots
# Plot for N = 20
df_log_gamma_eff_s20 <- as.data.frame(m_log_gamma_eff20_s$`true_value:m`)
levels(df_log_gamma_eff_s20$cond__) <- label_outdims
levels(df_log_gamma_eff_s20$effect2__) <- label_models

p_log_gamma_eff_s20 <- ggplot(data = df_log_gamma_eff_s20, aes(x = effect1__, y = estimate__, 
                                                           colour = effect2__, fill = effect2__)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_log_gamma_eff_s20$lower__, ymax = df_log_gamma_eff_s20$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = expression(log~gamma~(differenced)), colour = 'Models', fill = 'Models') +
  #guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#009E73", "#F0E442", "#D55E00")) + 
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#009E73", "#F0E442", "#D55E00")) #+ ggtitle('(a)')

# Plot for N = 50
df_log_gamma_eff_s50 <- as.data.frame(m_log_gamma_eff50_s$`true_value:m`)
levels(df_log_gamma_eff_s50$cond__) <- label_outdims
levels(df_log_gamma_eff_s50$effect2__) <- label_models[-1]

p_log_gamma_eff_s50 <- ggplot(data = df_log_gamma_eff_s50, aes(x = effect1__, y = estimate__, 
                                                             colour = effect2__, fill = effect2__)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_log_gamma_eff_s50$lower__, ymax = df_log_gamma_eff_s50$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = expression(log~gamma~(differenced)), colour = 'Models', fill = 'Models') +
  #guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6),
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_colour_manual(values = c("#E69F00", "#009E73", "#F0E442", "#D55E00")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#F0E442", "#D55E00")) #+ ggtitle('(b)')

# Plot for N = 20
df_log_gamma_eff_s200 <- as.data.frame(m_log_gamma_eff200_s$`true_value:m`)
levels(df_log_gamma_eff_s200$cond__) <- label_outdims
levels(df_log_gamma_eff_s200$effect2__) <- label_models[-1]

p_log_gamma_eff_s200 <- ggplot(data = df_log_gamma_eff_s200, aes(x = effect1__, y = estimate__, 
                                                             colour = effect2__, fill = effect2__)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_log_gamma_eff_s200$lower__, ymax = df_log_gamma_eff_s200$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = expression(log~gamma~(differenced)), colour = 'Models', fill = 'Models') +
  #guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6),
        axis.title.y = element_blank(), legend.position = 'none') +
  scale_colour_manual(values = c("#E69F00", "#009E73", "#F0E442", "#D55E00")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#F0E442", "#D55E00")) #+ ggtitle('(c)')

# Combine plots
p_log_gamma_eff <- ((p_log_gamma_eff20 + p_log_gamma_eff_s20) / 
                   (p_log_gamma_eff50 + p_log_gamma_eff_s50) / 
                   (p_log_gamma_eff200 + p_log_gamma_eff_s200)) + 
                   plot_layout(axis_titles = 'collect')
ggsave('se_log_gamma_eff.pdf',
       p_log_gamma_eff,
       dpi = 300,
       width = 40,
       height = 30,
       units = 'cm')
