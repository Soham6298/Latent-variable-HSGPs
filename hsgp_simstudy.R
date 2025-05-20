## Full simulation study for HSGPs and exact GPs. We generate data and fit models on cluster computers for 50 simulation trials.

#libraries
library(rstan)
library(posterior)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
# Import functions
source('hsgpfitfns.R')
# Compile stan models
hsgpmodel <- stan_model('hsgp_maternclass.stan')
gpmodel <- stan_model('gp_maternclass.stan')
# Set temporary folder (only to prevent cluster memory overflow)
#Sys.setenv(TMPDIR = "/mnt/volume")
#unlink(tempdir(), recursive = TRUE)
#tempdir(check = TRUE)

# Specify data generating conditions
N <- c(20, 50, 200)
dims <- c(5, 10, 20)
s_x <- 0.3
true_x_min <- 0
true_x_max <- 10
marginal_sd_params <- c(3, 0.25)
error_sd_params <- c(1, 0.05)
ls_params <- c(1, 0.25)
intc_params <- c(0, 10)
delta <- 1e-12
covfn <- 'se'

# Generate simulated data
n_sim <- 50
params <- list()
simdata <- list()
x_true <- list()
for (i in 1:n_sim) {
  set.seed(i)
  params[[i]] <- list()
  simdata[[i]] <- list() 
  x_true[[i]] <- list()
  for (j in 1:length(N)) {
    params[[i]][[j]] <- list()
    simdata[[i]][[j]] <- list()
    x_true[[i]][[j]] <- list()
    for (k in 1:length(dims)) {
      x_true[[i]][[j]][[k]] <- runif(N[j], true_x_min, true_x_max) 
      params[[i]][[j]][[k]] <- true_params_vary(n_obs = N[j], dims = dims[k],
                                                msd_pars = marginal_sd_params,
                                                esd_pars = error_sd_params,
                                                ls_pars = ls_params, 
                                                intc_pars = intc_params)
      simdata[[i]][[j]][[k]] <- gp_sim_data(n_obs = N[j], 
                                            dims = dims[k], 
                                            covfn = covfn,
                                            true_x = x_true[[i]][[j]][[k]], 
                                            s_x = s_x, 
                                            rho = params[[i]][[j]][[k]]$rho, 
                                            alpha = params[[i]][[j]][[k]]$alpha, 
                                            sigma = params[[i]][[j]][[k]]$sigma,
                                            intc = params[[i]][[j]][[k]]$intc)
    }
  }
}
# Save simulated data
saveRDS(params, 'hsgp_trueparams_se.rds')
saveRDS(simdata, 'hsgp_simdata_se.rds')

# Model specifications for fitting
x_min_model <- 0
x_max_model <- 10
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
marginal_sd_params_model <- c(3, 0.25) # alpha ~ N( , )
error_sd_params_model <- c(1, 0.25)    # sigma ~ N( , )
ls_params_model <- c(1, 0.05) # rho ~ N( , ) or InvGamma( , )
intc_params_model <- c(0, 10)
covfn_model <- 0  # 0 = SE; 1 = Matern3/2; 2 = Matern 5/2
# Set covariance function based scaler for HSGPs
if(covfn_model==1){
  const <- 3.42
}else if(covfn_model==2){
  const <- 2.65
}else {
  const <- 1.75
}
# Compute minimum number of basis functions
min_m <- ceiling(const * (5/4) * (x_max_model - x_min_model)/ls_params_model[1]) # mean(rho) = 0.3, S = 0.49 (half range of x), c = 5/4
# Add more basis function choices to test
M <- c(min_m, min_m + 4, min_m + 8)
# Set acceptance probability for MCMC 
adapt_delta_model <- 0.95

# Model fit and summary
## declare variable names
x_names <- list()
for (n in 1:length(N)) {
  x_names[[n]] <- sprintf('x[%s]', seq(1:N[n]))
}
rho_names <- list()
alpha_names <- list()
sigma_names <- list()
for (i in 1:length(dims)) {
  rho_names[[i]] <- sprintf('rho[%s]', seq(1:dims[i]))
  alpha_names[[i]] <- sprintf('alpha[%s]', seq(1:dims[i]))
  sigma_names[[i]] <- sprintf('sigma[%s]', seq(1:dims[i]))
}
# parallel for exact GP
out_gp <- list()
out_gp1 <- list()
cores = 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_exact = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  out_gp1[[i]] <- list()
  for (j in 1:length(dims)) {
    y <- simdata[[i]][[1]][[j]][,1:dims[j]]
    x_obs <- simdata[[i]][[1]][[j]]$x_obs
    x_true <- simdata[[i]][[1]][[j]]$x_true
    gp_fit <- gp_model(gpmodel = gpmodel, 
                       n_obs = N[1], 
                       dims = dims[j], 
                       outputs = y,
                       inputs = x_obs, 
                       latent_sd = s_x, 
                       latent_inputs = 1, 
                       covfn = covfn_model, 
                       is_vary = 1, 
                       is_corr = 0, 
                       rho_prior = rho_prior_model, 
                       ls_param = ls_params_model, 
                       x_min = x_min_model, 
                       x_max = x_max_model,
                       msd_param = marginal_sd_params_model, 
                       esd_param = error_sd_params_model, 
                       mean_intc = intc_params_model[1], 
                       sd_intc = intc_params_model[2], 
                       iter = 2000, 
                       warmup = 1000, 
                       chains = 1, 
                       cores = 1, 
                       init = 0, 
                       adapt_delta = adapt_delta_model)
    out_gp_x <- compare_summary(model = gp_fit, 
                                variable = x_names[[1]], 
                                true_variable = simdata[[i]][[1]][[j]]$x_true, 
                                n_obs = N[1], 
                                m_approx = 'exact', 
                                dims = dims[j], 
                                variable_class = 'x', 
                                model_name = 'gp', 
                                sim_id = i)
    out_gp_rho <- compare_summary(model = gp_fit, 
                                  variable = rho_names[[j]],
                                  true_variable = params[[i]][[1]][[j]]$rho, 
                                  n_obs = N[1], 
                                  m_approx = 'exact',
                                  dims = dims[j], 
                                  variable_class = 'rho', 
                                  model_name = 'gp', 
                                  sim_id = i)
    out_gp_alpha <- compare_summary(model = gp_fit, 
                                    variable = alpha_names[[j]],
                                    true_variable = params[[i]][[1]][[j]]$alpha,
                                    n_obs = N[1], 
                                    m_approx = 'exact',
                                    dims = dims[j], 
                                    variable_class = 'alpha', 
                                    model_name = 'gp',
                                    sim_id = i)
    out_gp_sigma <- compare_summary(model = gp_fit, 
                                    variable = sigma_names[[j]],
                                    true_variable = params[[i]][[1]][[j]]$sigma,
                                    n_obs = N[1],
                                    m_approx = 'exact',
                                    dims = dims[j], 
                                    variable_class = 'sigma', 
                                    model_name = 'gp',
                                    sim_id = i)
    out_gp1[[i]][[j]] <- rbind(out_gp_x, out_gp_rho, out_gp_alpha, out_gp_sigma)
  }
  out_gp[[i]] <- rbindlist(out_gp1[[i]])
}
compare_exact <- rbindlist(model_fit_exact)

# parallel for HSGP
out_hsgp <- list()
out_hsgp1 <- list()
out_hsgp2 <- list()
out_hsgp3 <- list()
cores = 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_hsgp = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  #set.seed(i)
  out_hsgp1[[i]] <- list()
  out_hsgp2[[i]] <- list()
  out_hsgp3[[i]] <- list()
  for (n in 1: length(N)) {
    out_hsgp2[[i]][[n]] <- list()
    out_hsgp3[[i]][[n]] <- list()
    for (j in 1:length(dims)) {
      y <- simdata[[i]][[n]][[j]][,1:dims[j]]
      x_obs <- simdata[[i]][[n]][[j]]$x_obs
      x_true <- simdata[[i]][[n]][[j]]$x_true
      out_hsgp3[[i]][[n]][[j]] <- list()
      for (k in 1:length(M)) {
        hsgp_fit <- hsgp_model(hsgpmodel = hsgpmodel, 
                               n_obs = N[n], 
                               m_obs = M[k], 
                               dims = dims[j], 
                               outputs = y, 
                               inputs = x_obs, 
                               latent_sd = s_x, 
                               latent_inputs = 1, 
                               c = 5/4, 
                               x_min = x_min_model, 
                               x_max = x_max_model,
                               covfn = covfn_model,
                               is_vary = 1, 
                               is_corr = 0, 
                               rho_prior = rho_prior_model, 
                               ls_param = ls_params_model, 
                               msd_param = marginal_sd_params_model, 
                               esd_param = error_sd_params_model, 
                               mean_intc = intc_params_model[1], 
                               sd_intc = intc_params_model[2], 
                               iter = 2000, 
                               warmup = 1000, 
                               chains = 1, 
                               cores = 1, 
                               init = 0, 
                               adapt_delta = adapt_delta_model)
        out_hsgp_x <- compare_summary(model = hsgp_fit, 
                                      variable = x_names[[n]], 
                                      true_variable = simdata[[i]][[n]][[j]]$x_true, 
                                      n_obs = N[n],
                                      m_approx = M[k], 
                                      dims = dims[j],
                                      variable_class = 'x', 
                                      model_name = 'hsgp', 
                                      sim_id = i)
        out_hsgp_rho <- compare_summary(model = hsgp_fit,
                                        variable = rho_names[[j]], 
                                        true_variable = params[[i]][[n]][[j]]$rho,
                                        n_obs = N[n], 
                                        m_approx = M[k], 
                                        dims = dims[j], 
                                        variable_class = 'rho', 
                                        model_name = 'hsgp',
                                        sim_id = i)
        out_hsgp_alpha <- compare_summary(model = hsgp_fit, 
                                          variable = alpha_names[[j]], 
                                          true_variable = params[[i]][[n]][[j]]$alpha, 
                                          n_obs = N[n], 
                                          m_approx = M[k], 
                                          dims = dims[j], 
                                          variable_class = 'alpha', 
                                          model_name = 'hsgp', 
                                          sim_id = i)
        out_hsgp_sigma <- compare_summary(model = hsgp_fit,
                                          variable = sigma_names[[j]], 
                                          true_variable = params[[i]][[n]][[j]]$sigma, 
                                          n_obs = N[n], 
                                          m_approx = M[k], 
                                          dims = dims[j], 
                                          variable_class = 'sigma',
                                          model_name = 'hsgp', 
                                          sim_id = i)
        out_hsgp3[[i]][[n]][[j]][[k]] <- rbind(out_hsgp_x, out_hsgp_rho, out_hsgp_alpha, out_hsgp_sigma)
      }
      out_hsgp2[[i]][[n]][[j]] <- rbindlist(out_hsgp3[[i]][[n]][[j]])
    }
    out_hsgp1[[i]][[n]] <- rbindlist(out_hsgp2[[i]][[n]])
  }
  out_hsgp[[i]] <- rbindlist(out_hsgp1[[i]])
}
compare_hsgp <- rbindlist(model_fit_hsgp)

## Bind to long table
compare_table <- rbind(compare_exact, compare_hsgp)

## comparison table
compare_table$data_id = paste0(compare_table$sim_id, '_', 
                               compare_table$m, '_',
                               compare_table$n, '_',
                               compare_table$d)
# Save simulation study results
saveRDS(compare_table, 'hsgp_simout_se.rds')
