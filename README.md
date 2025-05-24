# Hilbert space methods for approximating multi-output latent variable Gaussian processes
This is a repository for the codes used in the paper: Hilbert space methods for approximating multi-output latent variable Gaussian processes (https://doi.org/10.48550/arXiv.2505.16919). The codes are arranged as follows. Check 'hsgpfitfns.R' to call necessary functions that can be used in all of the R scripts (for e.g: to customise any simulation conditions, fit models, prepare summary plots).

## Models
The approximate multi-output latent variable Hilbert space Gaussian process (HSGPs) for the Matern class of covariance functions model develpment are given as 'hsgp_maternclass.stan'. The exact Gaussian process (GP) counterparts are provided in 'gp_maternclass.stan'. The variational latent GP (VIGP) based on pyro-GPLVM  module (https://pyro.ai/examples/gplvm.html) is provided in the 'vigp_simstudy.py' and 'vigp_casestudy.py' along with their respective studies.

## Simulation Study
See 'hsgp_simstudy.R' for the entire simulation study design. The script provides an output dataframe with all the results for post-processing. Export the generated data from the R script using 'exportSimdata_R2P_cluster.R' before running experiments for VIGP using 'vigp_simstudy.py'. Our simulation experiment results are provided in the 'hsgp and exact gp sim results' and 'vigp' sim results' folders for the respective models that can be post-analysed using the appropriate R scripts.

## Model convergence
Check model convergence for all our simulation experiments using the 'ModelConvergence.R' script.

## Model Calibration
We use SBC (https://projecteuclid.org/journals/bayesian-analysis/volume-20/issue-2/Simulation-Based-Calibration-Checking-for-Bayesian-Computation--The-Choice/10.1214/23-BA1404.full) to check for our model calibration. See 'sbc_ecdf_plots.R' and 'log_gamma_summary.R' to generate rank ECDF graphical tests and log gamma scores respectively for model calibration using the simulation results.

## Latent variable estimation
Check latent variable estimation accuracy for various simulation studies provided (or create your own) and compare between HSGPs, exact GPs and VIGPs using the 'LatentxRecovery.R' script.

## GP hyperparameter estimates
Check the GP hyperparameter estimates and recovery of true simulated values using 'HyperparameterRecovery.R'.

## Case study
We apply HSGPs on a single-cell biology data depicting cyclic cellular phases (https://www.nature.com/articles/s41586-021-03232-9) to esimate latent cellular ordering. We additionally compare the results with VIGPs (fitted using 'vigp_casestudy.py'). See 'CaseStudy.R' to re-analyse HSGPs and the comparative study using the data provided in 'case study data' folder. For more details, check the 'Real-World Case Study' section of our paper.
