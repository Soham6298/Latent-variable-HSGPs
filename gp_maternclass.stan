functions {
  // base covariance fn
  
  // SE
    matrix se(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = exp(r * square(x[i] - x[j])) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  // Matern 3/2
  matrix m32(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  //Matern 5/2 
   matrix m52(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                    ((5 * square(x[i] - x[j]))/ (3 * sq_rho))) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
}

data {
  // Sample size
  int<lower=1> N;
  // Number of output dimensions
  int<lower=1> D;
  // Output
  matrix[N, D] y; 
  // Prior measurement SD (assumed to be known)
  real<lower=0> s;
  // Prior for latent x
  vector[N] inputs;
  real x_min;      // lower bound for true_x prior
	real x_max;      // upper bound for true_x prior
  // Model conditions (specify which modifications to use ranging from standard GPs to DGP-LVM)
  int<lower=0, upper=1> is_vary; // 0 = constant rho and alpha for each dims; 1 = varying params
  int<lower=0, upper=1> is_corr; // 0 = no correlation b/w dims; 1 = correlated dims
  int<lower=0, upper=2> covfn; //0 = SE; 1 = M3/2; 2 = M5/2
  int<lower=0, upper=1> latent; //0 = obs inputs; 1 = latent inputs
  int<lower=0, upper=1> rho_prior; //0 = normal; 1 = invgamma;
  // prior hyperparams specification (only for two parameter dist families)
  vector[2] ls_param; 
  vector[2] msd_param;
  vector[2] esd_param;
  // Input data specific mean and sd for the y and y'
  real mean_intc;
  real sd_intc;
}

transformed data {
// add constant to the diagonal of the covariance matrix for computational stability
  real delta = 1e-6;
}

parameters {
  // GP Length scale parameter (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] rho_temp;
  // GP marginal SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] alpha_temp;
  // Gaussian link parameter
  matrix[N, D] eta;
  // Between dimension correlation parameter (temporary storage)
  cholesky_factor_corr[is_corr==1 ? D:0] L_omega_temp;
  // Error SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] sigma_temp;
  // only for latent variable models
  vector<lower=x_min,upper=x_max>[latent==1 ? N:0] x_temp;
  array[D] real intercept;
}

transformed parameters {
  // Latent x
  // add bounds on x if known for theoretical reasons (e.g. lower=0)
  vector[N] x;
  if (latent == 1) {
	  x = x_temp;
	} else{
	  x = inputs;
	}
  // GP f()
  matrix[N, D] f;
  // Model condition adjustments (temporary to actual parameters)
  vector<lower=0>[D] rho;
  vector<lower=0>[D] alpha;
  vector<lower=0>[D] sigma;
  cholesky_factor_corr[D] L_omega;
  // if params are constant for dims, it will be repeated
  if (is_vary == 1) {
    rho = rho_temp;
    alpha = alpha_temp;
    sigma = sigma_temp;
  } else {
    for (k in 1:D) {
      rho[k] = rho_temp[1];
      alpha[k] = alpha_temp[1];
      sigma[k] = sigma_temp[1];
    }
  }
  // if no correlation b/w dims, corr mat will be replaced by identity matrix
  if (is_corr==0) {
    L_omega = diag_matrix(rep_vector(1, D));
  } else {
    L_omega = L_omega_temp;
  }
  // Computing covariance matrix for standard GP
  if (covfn == 0) {
    for (k in 1:D) {
      matrix[N, N] K = se(x, alpha[k], rho[k], delta);
      f[, k] = K * eta[, k];
    }
  } else if (covfn == 1) {
    for (k in 1:D) {
      matrix[N, N] K = m32(x, alpha[k], rho[k], delta);
      f[, k] = K * eta[, k];
    } 
  } else if (covfn == 2) {
    for (k in 1:D) {
      matrix[N, N] K = m52(x, alpha[k], rho[k], delta);
      f[, k] = K * eta[, k];
    } 
  }
    // For correlated outputs
    f = f * L_omega';
}


model {
  // Priors
  if (rho_prior == 0) {
    rho_temp ~ normal(ls_param[1], ls_param[2]); 
  } else if (rho_prior == 1) {
    rho_temp ~ inv_gamma(ls_param[1], ls_param[2]);
  }
  alpha_temp ~ normal(msd_param[1], msd_param[2]);
  sigma_temp ~ normal(esd_param[1], esd_param[2]);
  L_omega_temp ~ lkj_corr_cholesky(1);
  for (k in 1:D) {
    to_vector(eta[, k]) ~ std_normal();
  }
  for (k in 1:D) {
    intercept[k] ~ normal(mean_intc, sd_intc);
  }
  // Likelihood
  if (latent) {
    inputs ~ normal(x_temp, s);
	}
	// set prior on x to match data generating process 
	x_temp ~ uniform(x_min, x_max);
  for (k in 1:D) {
    y[, k] ~ normal(intercept[k] + f[, k], sigma[k]);
  }
}
