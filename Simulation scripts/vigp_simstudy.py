# Import dependencies
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy as scp
from joblib import Parallel, delayed
import multiprocessing
import time

import torch
import pathlib
from torch.nn import Parameter

import pyro
import pyro.contrib.gp as gp
import pyro.distributions as dist
import pyro.ops.stats as stats
from pyro.nn import PyroSample
from torch.nn import Parameter

smoke_test = ('CI' in os.environ)  # ignore; used to check code integrity in the Pyro repo
assert pyro.__version__.startswith('1.9.1')
pyro.set_rng_seed(123)

# First create the simulation data from the hsgp_simstudy.R script
# Then export the data using the exportSimdata_R2P_cluster.R script

# Ensure simdata structure is correct and properly reset
simdata = []

# Read all simdata.csv from path for each N and dim
N = [20, 50, 200]
dims = [5, 10, 20]
for n in N:
    simdata_dim = []  # Reset for each sample size
    for dim in dims:
        simdatapath = glob.glob(f'Simdata{n}N/Simdata{dim}D/*.csv')
        simdatapath = sorted(simdatapath)
        data = [pd.read_csv(file, index_col=0) for file in simdatapath]
        simdata_dim.append(data)  # Store as list of DataFrames
    simdata.append(simdata_dim)

simdata[0][0][0].head()

## Model fit function
def gplvmfit(output, input_prior):
    pyro.clear_param_store()
    # Convert input data to tensors
    data = torch.tensor(output.values, dtype=torch.get_default_dtype())
    y = data.t()  # Transpose to shape [num_data_points, num_features]
    
    time = torch.tensor(input_prior.values.squeeze(), dtype=torch.get_default_dtype())
    
    # Dynamically update the size of X_prior_mean based on output and input_prior lengths
    X_prior_mean = torch.zeros(y.size(1), 1, dtype=torch.get_default_dtype()) 
    X_prior_mean[:, 0] = time[:y.size(1)]  # Ensure that X_prior_mean matches the size of y
    
    # Define kernel
    kernel = gp.kernels.RBF(input_dim=1)
    #kernel = gp.kernels.Matern32(input_dim=1)
    #kernel = gp.kernels.Matern52(input_dim=1)
    
    # Use Parameter so PyTorch can track it for optimization
    X = Parameter(X_prior_mean.clone(), requires_grad=True)
    
    # Sampling inducing points using random resampling
    Xu = X_prior_mean[torch.randperm(X_prior_mean.size(0))[:10]]
    
    # Create the Sparse GP model (lower jitter will probably cause positive-definiteness error)
    gplvm = gp.models.SparseGPRegression(X, y, kernel, Xu, jitter=1e-3)
    
    # Pyro sample for X with a normal prior
    gplvm.X = pyro.nn.PyroSample(dist.Normal(X_prior_mean, 0.3).to_event())
    ## Setting custom priors raise warnings
    #gplvm.kernel.lengthscale = PyroSample(dist.HalfNormal(1, 0.05))
    #gplvm.kernel.variance = PyroSample(dist.HalfNormal(3, 0.25))
    #gplvm.noise = PyroSample(dist.HalfNormal(1, 0.25))
    
    # Automatic guide for variational inference
    gplvm.autoguide("X", dist.Normal)
    
    # Train the model
    losses = gp.util.train(gplvm, num_steps=4000)
    
    return gplvm

# Sim trial function
def simtrial(i, n, j, dim):
    cols = [f'X{k}' for k in range(1, dim + 1)]

    # Extract the relevant data from simdata
    X = simdata[n][j][i].loc[:, cols]
    y = simdata[n][j][i][['x_obs']]

    # Debugging: print the shapes of X and y to ensure they match where expected
    print(f"Shape of X (features): {X.shape}")
    print(f"Shape of y (observations): {y.shape}")

    # Ensure that the number of rows (observations) match
    if X.shape[0] != y.shape[0]:
        raise ValueError(f"Mismatch in number of rows: X has {X.shape[0]} rows, y has {y.shape[0]} rows")
    fitmodel = gplvmfit(X, y)
    #fitmodel.autoguide("X", dist.Normal)    
    return fitmodel

# Fit models in parallel
n_sim = 50
num_cores = 50

# Main list to store fitted models
fittedmodels = []

# Loop through each sample size (N values)
for n, obs in enumerate(N):
    print(f"Processing for {obs} sample size...")
    fittedmodels_dim = []  # Reset for each sample size
    
    # Loop through each dimension for the given sample size
    for j, dim in enumerate(dims):
        print(f"Processing for {dim} dimensions...")

        # Perform parallel simulation trials
        fit = Parallel(n_jobs=num_cores)(delayed(simtrial)(i, n, j, dim) for i in range(n_sim))
        # Append the fits for this dimension
        fittedmodels_dim.append(fit)
        
        # Verify the number of fits
        print(f"Finished processing for {dim}D. Number of fits: {len(fit)}")
    
    # Append the data for all dimensions for the current sample size
    fittedmodels.append(fittedmodels_dim)
    print(f"Finished processing for {obs}N. Number of dimension sets: {len(fittedmodels_dim)}")

# Final check: print total counts
print(f"Total number of sample size sets: {len(fittedmodels)}")

# Arrange summary outputs
X_mean, X_sd, Rho, Alpha, Sigma, X_true = [], [], [], [], [], []

for n in range(len(N)):
    X_mean_n, X_sd_n, Rho_n, Alpha_n, Sigma_n, X_true_n = [], [], [], [], [], []
    
    for j in range(len(dims)):
        X_mean_j, X_sd_j, Rho_j, Alpha_j, Sigma_j, X_true_j = [], [], [], [], [], []
        
        for i in range(n_sim):
            # Ensure mode setting is correct
            fittedmodels[n][n][i].mode = "guide"
        
            # Extract and convert the required attributes
            X_mean_temp = fittedmodels[n][j][i].X_loc.detach().cpu().numpy()  # Ensure data is on CPU
            X_sd_temp = fittedmodels[n][j][i].X_scale.detach().cpu().numpy()
            Rho_temp = fittedmodels[n][j][i].kernel.lengthscale.detach().cpu().numpy()
            Alpha_temp = fittedmodels[n][j][i].kernel.variance.detach().cpu().numpy()
            Sigma_temp = fittedmodels[n][j][i].noise.detach().cpu().numpy()
            X_true_temp = np.array(simdata[n][j][i][['x_true']].values.squeeze())
        
            # Append temporary values to lists for current n
            X_mean_j.append(X_mean_temp)
            X_sd_j.append(X_sd_temp)
            Rho_j.append(Rho_temp)
            Alpha_j.append(Alpha_temp)
            Sigma_j.append(Sigma_temp)
            X_true_j.append(X_true_temp)
    
        # Append the lists for the current n
        X_mean_n.append(X_mean_j)
        X_sd_n.append(X_sd_j)
        Rho_n.append(Rho_j)
        Alpha_n.append(Alpha_j)
        Sigma_n.append(Sigma_j)
        X_true_n.append(X_true_j)
    # Append the lists for the current n
    X_mean.append(X_mean_n)
    X_sd.append(X_sd_n)
    Rho.append(Rho_n)
    Alpha.append(Alpha_n)
    Sigma.append(Sigma_n)
    X_true.append(X_true_n)

# Define correct N values for reference
N_values = [20, 50, 200]  # Sample sizes
Dims_values = [5, 10, 20]   # Output dimensions
Sim = 50  # Number of trials

# Initialize lists to hold flattened data
X_mean_flat, X_sd_flat, X_true_flat = [], [], []
N_labels, D_labels, Sim_labels, Obs_labels = [], [], [], []

# Flatten the nested lists while keeping track of i, j, k indices
for n_idx, (N, mean_n, sd_n, true_n) in enumerate(zip(N_values, X_mean, X_sd, X_true)):  
    for j_idx, (D, mean_dim, sd_dim, true_dim) in enumerate(zip(Dims_values, mean_n, sd_n, true_n)):
        for sim_idx in range(Sim):  # 50 simulations per (N, D)

            # Check available observations (should match N)
            num_obs = len(mean_dim[sim_idx])  
            if num_obs != N:
                print(f"⚠️ Warning: N={N}, Sim={sim_idx}, Expected Obs={N}, Found Obs={num_obs}")

            for obs_idx in range(min(N, num_obs)):  # Ensure Obs goes from 1 to N
                X_mean_flat.append(mean_dim[sim_idx][obs_idx])
                X_sd_flat.append(sd_dim[sim_idx][obs_idx])
                X_true_flat.append(true_dim[sim_idx][obs_idx])

                # Append correct labels
                N_labels.append(N)  # Correct Sample Size
                D_labels.append(D)  # Output Dimensions
                Sim_labels.append(sim_idx + 1)  # Simulation number
                Obs_labels.append(obs_idx + 1)  # 1-based Obs index

# Create DataFrame with correct labels
df = pd.DataFrame({
    'N': N_labels,
    'Dims': D_labels,
    'Sim': Sim_labels,
    'Obs': Obs_labels,  # Ensures Obs matches 1 → N correctly
    'X_true': X_true_flat,
    'X_mean': X_mean_flat,
    'X_sd': X_sd_flat
})
# View the DataFrame
print(df.head())  # Printing the first few rows of the DataFrame
# Save to CSV
df.to_csv('latentXvalues.csv', index=False)
# Read the CSV file and clean up any brackets [] in the data
df_cleaned = pd.read_csv("latentXvalues.csv").apply(lambda x: x.replace(r"[\[\]]", "", regex=True))

# Save the cleaned DataFrame to a new CSV file
df_cleaned.to_csv("GPLVMpyOut_se.csv", index=False)