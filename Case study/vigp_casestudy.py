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
pyro.set_rng_seed(1)

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
    
    # Create the Sparse GP model
    gplvm = gp.models.SparseGPRegression(X, y, kernel, Xu, jitter=1e-3)
    
    # Pyro sample for X with a normal prior
    gplvm.X = pyro.nn.PyroSample(dist.Normal(X_prior_mean, 0.03).to_event())
    
    # Automatic guide for variational inference
    gplvm.autoguide("X", dist.Normal)
    
    # Train the model
    losses = gp.util.train(gplvm, num_steps=4000)
    
    return gplvm

# import data from case study data folder
data = pd.read_csv('case_study_hsgp_data.csv', index_col=0)
data.head()
data.columns = data.columns.str.strip()
print(data.columns)

x = data[['t']]  # if 't' is correct after checking
y = data.drop(columns=['t', 'p'])
fitmodel = gplvmfit(y, x)

x_mean_post = fitmodel.X_loc.detach().cpu().numpy().flatten()
x_sd_post = fitmodel.X_scale.detach().cpu().numpy().flatten()
t = data['t'].values
print(len(t), len(x_mean_post), len(x_sd_post))

df = pd.DataFrame({
    't': t,
    'X_mean': x_mean_post,
    'X_sd': x_sd_post
})
# View the DataFrame
print(df.head())  # Printing the first few rows of the DataFrame
# Save to CSV
df.to_csv('latentXvalues.csv', index=False)
# Read the CSV file and clean up any brackets [] in the data
df_cleaned = pd.read_csv("latentXvalues.csv").apply(lambda x: x.replace(r"[\[\]]", "", regex=True))

# Save the cleaned DataFrame to a new CSV file
df_cleaned.to_csv("GPLVMpyOut_se_casestudy.csv", index=False)