{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "daa26fef-5ead-44ea-b963-03c37f2e4d91",
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\n",
    "import glob\n",
    "import matplotlib.pyplot as plt\n",
    "import pandas as pd\n",
    "import numpy as np\n",
    "import scipy as scp\n",
    "from joblib import Parallel, delayed\n",
    "import multiprocessing\n",
    "import time\n",
    "\n",
    "import torch\n",
    "import pathlib\n",
    "from torch.nn import Parameter\n",
    "\n",
    "import pyro\n",
    "import pyro.contrib.gp as gp\n",
    "import pyro.distributions as dist\n",
    "import pyro.ops.stats as stats\n",
    "from pyro.nn import PyroSample\n",
    "from torch.nn import Parameter\n",
    "\n",
    "smoke_test = ('CI' in os.environ)  # ignore; used to check code integrity in the Pyro repo\n",
    "assert pyro.__version__.startswith('1.9.1')\n",
    "pyro.set_rng_seed(1)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c2a70490-1398-4228-8452-2ca2b7216f42",
   "metadata": {},
   "outputs": [],
   "source": [
    "## Model fit function\n",
    "def gplvmfit(output, input_prior):\n",
    "    pyro.clear_param_store()\n",
    "    # Convert input data to tensors\n",
    "    data = torch.tensor(output.values, dtype=torch.get_default_dtype())\n",
    "    y = data.t()  # Transpose to shape [num_data_points, num_features]\n",
    "    \n",
    "    time = torch.tensor(input_prior.values.squeeze(), dtype=torch.get_default_dtype())\n",
    "    \n",
    "    # Dynamically update the size of X_prior_mean based on output and input_prior lengths\n",
    "    X_prior_mean = torch.zeros(y.size(1), 1, dtype=torch.get_default_dtype()) \n",
    "    X_prior_mean[:, 0] = time[:y.size(1)]  # Ensure that X_prior_mean matches the size of y\n",
    "    \n",
    "    # Define kernel \n",
    "    kernel = gp.kernels.RBF(input_dim=1)\n",
    "    #kernel = gp.kernels.Matern32(input_dim=1)\n",
    "    #kernel = gp.kernels.Matern52(input_dim=1)\n",
    "    \n",
    "    # Use Parameter so PyTorch can track it for optimization\n",
    "    X = Parameter(X_prior_mean.clone(), requires_grad=True)\n",
    "    \n",
    "    # Sampling inducing points using random resampling\n",
    "    Xu = X_prior_mean[torch.randperm(X_prior_mean.size(0))[:10]]\n",
    "    \n",
    "    # Create the Sparse GP model\n",
    "    gplvm = gp.models.SparseGPRegression(X, y, kernel, Xu, jitter=1e-3)\n",
    "    \n",
    "    # Pyro sample for X with a normal prior\n",
    "    gplvm.X = pyro.nn.PyroSample(dist.Normal(X_prior_mean, 0.03).to_event())\n",
    "    \n",
    "    # Automatic guide for variational inference\n",
    "    gplvm.autoguide(\"X\", dist.Normal)\n",
    "    \n",
    "    # Train the model\n",
    "    losses = gp.util.train(gplvm, num_steps=4000)\n",
    "    \n",
    "    return gplvm"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2a094f6c-51f9-4b7d-862e-499ae4ad78aa",
   "metadata": {},
   "outputs": [],
   "source": [
    "# import data from case study data folder\n",
    "data = pd.read_csv('case_study_hsgp_data.csv', index_col=0)\n",
    "data.head()\n",
    "data.columns = data.columns.str.strip()\n",
    "print(data.columns)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "bfe39f23-b063-4c1a-80df-e66a1d393984",
   "metadata": {},
   "outputs": [],
   "source": [
    "x = data[['t']]  # if 't' is correct after checking\n",
    "y = data.drop(columns=['t', 'p'])\n",
    "fitmodel = gplvmfit(y, x)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "e3edb9c0-0c5e-4b5a-8fac-0486bb4e72be",
   "metadata": {},
   "outputs": [],
   "source": [
    "x_mean_post = fitmodel.X_loc.detach().cpu().numpy().flatten()\n",
    "x_sd_post = fitmodel.X_scale.detach().cpu().numpy().flatten()\n",
    "t = data['t'].values\n",
    "print(len(t), len(x_mean_post), len(x_sd_post))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7050bdca-1544-413a-9efa-dee16b9c4382",
   "metadata": {},
   "outputs": [],
   "source": [
    "df = pd.DataFrame({\n",
    "    't': t,\n",
    "    'X_mean': x_mean_post,\n",
    "    'X_sd': x_sd_post\n",
    "})\n",
    "# View the DataFrame\n",
    "print(df.head())  # Printing the first few rows of the DataFrame\n",
    "# Save to CSV\n",
    "df.to_csv('latentXvalues.csv', index=False)\n",
    "# Read the CSV file and clean up any brackets [] in the data\n",
    "df_cleaned = pd.read_csv(\"latentXvalues.csv\").apply(lambda x: x.replace(r\"[\\[\\]]\", \"\", regex=True))\n",
    "\n",
    "# Save the cleaned DataFrame to a new CSV file\n",
    "df_cleaned.to_csv(\"GPLVMpyOut_se_casestudy.csv\", index=False)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
