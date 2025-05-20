## We estimate the perform model evaluation for exact and HSGPs based on GP hyperparameter recovery from simulation studies
## Need to load the output dataframe from the simulation study code before generating figures here
## We only show hyperparameter recovery for the most complicated case of N = 200 and D = 20. The other cases are qualitatively similar

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
source('hsgpfitfns.R')
# Summary model
## load simulation study output
compare_table <- readRDS('gprout/hsgp_simout_se.rds')
# Filter for N = 200; D = 20
compare_n200 <- subset(compare_table, model_name=='hsgp' & n=='200' & d=='20')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
# Designate levels properly to set figure legends and guides
levels(compare_table$m)
compare_table <- compare_table %>% 
  mutate(m = recode(m, "exact" = "Exact", "34" = "HSGP(34)", "38" = "HSGP(38)", "42" = "HSGP(42)"))
compare_table <- compare_table %>% 
  mutate(d = recode(d, "5" = 'D = 5', '10' = 'D = 10', '20' = 'D = 20'))
compare_table <- compare_table %>% 
  mutate(n = recode(n, "20" = 'N = 20', '50' = 'N = 50', '200' = 'N = 200'))
# Set data for generating figure
compare_table$data_id <- as.factor(compare_table$data_id)
compare_pars <- subset(compare_table, class != 'x')
compare_pars <- subset(compare_pars, d == 'D = 20')
# Generate figure
p_pars_rmse <- ggplot(compare_pars, aes(x = class, y = rmse, colour = m)) +
  theme_bw(base_size=20,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  facet_wrap(~n) +
  labs(x = 'Parameters', y = 'RMSE', colour = 'Models') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#009E73", "#F0E442", "#D55E00" )) + ggtitle('')

ggsave('m52_hyperparams.pdf',
       p_pars_rmse,
       dpi = 300,
       width = 30,
       height = 10,
       units = 'cm')
