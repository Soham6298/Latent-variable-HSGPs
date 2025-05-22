# Hilbert space methods for approximating multi-output latent variable Gaussian processes
This is a repository for the codes used in the paper: Hilbert space methods for approximating multi-output latent variable Gaussian processes. The codes are arranged as follows. Check 'hsgpfitfns.R' to call necessary functions that can be used in all of the R scripts (for e.g: to customise any simulation conditions, fit models, prepare summary plots).

## Models
The approximate multi-output latent variable HSGPs for the Matern class of covariance functions model develpment are given as 'hsgp_maternclass.stan'. The exact GP counterparts are provided in 'gp_maternclass.stan'. The variational latent GP based on pyro-GPLVM  module (https://pyro.ai/examples/gplvm.html) is provided in the 'vigp_simstudy.py' and 'vigp_casestudy.py' along with their respective studies.

## Simulation Study
See 'hsgp_simstudy.R' for the entire simulation study design. The script provides an output dataframe with all the results for post-processing. The simulation experiments results are provided in the 'hsgp and exact gp sim results' and 'vigp' sim results' folders for the respective models.

## Model convergence
Check model convergence for all our simulation experiments using the 'ModelConvergence.R' script.

## Model Calibration
We use SBC (https://projecteuclid.org/journals/bayesian-analysis/volume-20/issue-2/Simulation-Based-Calibration-Checking-for-Bayesian-Computation--The-Choice/10.1214/23-BA1404.full) to check for our model calibration
