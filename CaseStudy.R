# We showcase HSGP for a single-cell RNA sequencing data
# The data cellcycle is obtained from Cytopath (https://doi.org/10.1016/j.crmeth.2022.100359)
# We provide recovery of latent input (pseudotime estimation) for HSGP and Pyro-VI with Squared Exponential covariance functions

#libraries
library(rstan)
library(bayesplot)
library(ggplot2)
library(posterior)
library(dplyr)
library(brms)
library(grid)
library(gridExtra)
library(patchwork)
# Source functions
source('hsgpfitfns.R')
stanmodel <- stan_model('hsgp_maternclass.stan')
# Read data
cellcycle <- readRDS('case study data/DGPLVMcasestudydata.Rdata')
cellcycle_obs <- read.csv('case study data/obs.csv')
cellcycledata <- data.frame(cellcycle, cell_cycle_phase = cellcycle_obs$cell_cycle_phase)
# Data cleanup
cellcycledata <- na.omit(cellcycledata)
# data for model
y <- cellcycledata[,3:14]
t <- cellcycledata$cell_cycle_hrs
p <- cellcycledata$cell_cycle_phase
# prepare final data for pyro-VI fitting
finaldata <- data.frame(y, t, p)
write.csv(finaldata, 'case study data/case_study_hsgp_data.csv')
# set for prior measurement SD for latent x
s_x <- 0.03
# Number of cells (sample size)
N <- nrow(y)
# Number of genes (output dimensions)
D <- ncol(y)
# model fitting
mean_obs <- rep(NA, D)
sd_obs <- rep(NA, D)
for(j in 1:D) {
    mean_obs[j] <- mean(y[1:N,j])
    sd_obs[j] <- sd(y[1:N,j])
}
# Model specifications
true_x_max <- 1
true_x_min <- 0
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
marginal_sd_params_model <- c(14, 3.5)  # alpha ~ N( , )
error_sd_params_model <- c(7, 3.5)     # sigma ~ N( , )
ls_params_model_se <- c(0.4, 0.1) # rho ~ N( , ) or InvGamma( , )
min_m <- ceiling(1.75 * (5/4) * (true_x_max - true_x_min)/(ls_params_model_se[1]))
intc_prior_model <- c(mean(mean_obs), mean(sd_obs))
latent_model <- 1 # 1 for latent inputs, 0 for manifest
adapt_delta_model <- 0.95
# Model fit
# SE
hsgp_se_fit <- hsgp_model(stanmodel,
                   n_obs = N,
                   m_obs = min_m,  
                   x_min = true_x_min, 
                   x_max = true_x_max, 
                   c = 5/4, 
                   dims = D, 
                   outputs = y, 
                   inputs = t, 
                   latent_sd = s_x, 
                   latent_inputs = latent_model, 
                   covfn = 0, #covfn_model <- 0  # 0 - SE; 1 = Matern3/2; 2 = Matern 5/2
                   rho_prior = rho_prior_model,
                   ls_param = ls_params_model_se, 
                   msd_param = marginal_sd_params_model, 
                   esd_param = error_sd_params_model, 
                   mean_intc = intc_prior_model[1],
                   sd_intc = intc_prior_model[2], 
                   is_vary = 1, 
                   is_corr = 1, 
                   iter = 2000, 
                   warmup = 1000, 
                   chains = 2, 
                   cores = 2, 
                   init = 0, 
                   adapt_delta = adapt_delta_model)

print(hsgp_se_fit, pars = c('rho','alpha', 'sigma', 'x'))
# HSGP vs cellhours figure
hsgp_se_draws <- as_draws_df(hsgp_se_fit)
x_names <- sprintf('x[%s]', seq(1:N))
hsgp_se_draws_x <- subset_draws(hsgp_se_draws, variable = x_names, chain = 1)
hsgp_se_summary <- summarise_draws(hsgp_se_draws_x)
hsgp_se_summary$obs_t <- t 
hsgp_se_summary$phases <- p
hsgp_se_summary$label <- 'HSGP vs Cell hours'
hsgp_se_plot <- ggplot(hsgp_se_summary, aes(x = obs_t, y = mean - obs_t, colour = phases)) + 
  theme_bw(base_size = 25, base_family = 'Times') +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  ylim(c(-0.35, 0.35)) +
  facet_wrap(~label) +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#009E73")) +
  labs(x = 'Cell hours', y = 'Difference', colour = 'Phases')

# pyro-VI vs cellhours figure
pyro_casestudy <- read.csv('gplvmout/GPLVMpyOut_se_casestudy.csv') ## output from fitted pyro-VI
pyro_casestudy$phases <- p
pyro_casestudy$label <- 'PyroVI vs Cell hours'
pyro_se_plot <- ggplot(pyro_casestudy, aes(x = t, y = X_mean - t, colour = phases)) + 
  theme_bw(base_size = 25, base_family = 'Times') +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  ylim(c(-0.35, 0.35)) +
  facet_wrap(~label) +
  labs(x = 'Cell hours', y = 'Difference', colour = 'Phases') +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#009E73")) 
# pyro-VI vs HSGP figure
pyro_hsgp_diff <- data.frame(obs_t = hsgp_se_summary$obs_t, 
                             phases = hsgp_se_summary$phases, 
                             hsgp_pyro_diff = hsgp_se_summary$mean - pyro_casestudy$X_mean, 
                             label = 'HSGP vs PyroVI')
pyro_hsgp_se_plot <- ggplot(pyro_hsgp_diff, aes(x = obs_t, y = hsgp_pyro_diff, colour = phases)) + 
  theme_bw(base_size = 25, base_family = 'Times') +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  ylim(c(-0.35, 0.35)) +
  facet_wrap(~label) +
  labs(x = 'Cell hours', y = 'Difference', colour = 'Phases') +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#009E73")) 
# combine all sub-figures
case_study_plot <- hsgp_se_plot + pyro_se_plot + pyro_hsgp_se_plot + plot_layout(axes = 'collect', guides = 'collect') &
                      theme(legend.position = 'bottom')
ggsave('hsgp_case_study_plot.pdf',
       case_study_plot,
       dpi = 300,
       width = 45,
       height = 15,
       units = 'cm')

## Figures for GP hyperparameter estimation
hsgp_se_draws_pars <- subset_draws(hsgp_se_draws, variable = c('rho', 'alpha', 'sigma'), chain = 1)
hsgp_se_summary_pars <- summarise_draws(hsgp_se_draws_pars)
var_name <- c(rep('rho', D), rep('alpha', D), rep('sigma', D))
output_dim <- rep(seq(1:D), 3) # Based on number of hyperparameters
# Design the dataframe for plot
hsgp_se_summary_pars$var_name <- var_name
hsgp_se_summary_pars$output_dim <- as.factor(output_dim)
hsgp_se_summary_rho <- hsgp_se_summary_pars[var_name=='rho',]
hsgp_se_summary_rho$ftitle <- 'rho'
hsgp_se_summary_alpha <- hsgp_se_summary_pars[var_name=='alpha',]
hsgp_se_summary_alpha$ftitle <- 'alpha'
hsgp_se_summary_sigma <- hsgp_se_summary_pars[var_name=='sigma',]
hsgp_se_summary_sigma$ftitle <- 'sigma'

# GP length-scale figure
hsgp_se_rho_plot <- ggplot(hsgp_se_summary_rho, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 40, base_family = 'Times') +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1) + 
  labs(x = '', y = '') +
  theme(axis.title = element_blank()) +
  facet_wrap(~ftitle) +
  scale_y_log10()
  
# GP marginal SD figure
hsgp_se_alpha_plot <- ggplot(hsgp_se_summary_alpha, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 40, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = '', y = '') + 
  facet_wrap(~ftitle) +
  scale_y_log10()

# Error SD figure
hsgp_se_sigma_plot <- ggplot(hsgp_se_summary_sigma, aes(x = output_dim, y = mean)) +   
  theme_bw(base_size = 40, base_family = 'Times') +
  geom_point(size = 3, position = position_dodge(0.7)) +
  geom_errorbar(aes(ymin = q5, ymax = q95), linewidth = 1, position = position_dodge(0.7)) + 
  labs(x = '', y = '') + 
  facet_wrap(~ftitle) +
  scale_y_log10()


# Combine al sub-figures
hsgp_cellcycle_pars <- (hsgp_se_rho_plot + hsgp_se_alpha_plot + hsgp_se_sigma_plot) +
  plot_layout(axis_titles = 'collect', guides = 'collect') &
  theme(axis.text = element_text(size = 40))
hsgp_cellcycle_pars
gt <- patchwork::patchworkGrob(hsgp_cellcycle_pars)
g <- arrangeGrob(gt, left = textGrob("Estimate", rot = 90, gp = gpar(fontsize=40, fontfamily='Times')), 
                  bottom = textGrob("Gene", vjust = -0.2, gp = gpar(fontsize=40, fontfamily='Times')))
ggsave('hsgp_case_study_pars.pdf',
       g,
       dpi = 300,
       width = 80,
       height = 20,
       units = 'cm')
