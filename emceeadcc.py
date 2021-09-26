import numpy as np
from theano import tensor as T
import matplotlib.pyplot as plt
from deconv.imports import load_tables, load_figures
import pymc3 as pm
import arviz as az

def makeFigure():
    M = pm.Model()

    A_antiD, _, _, _ = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    mean_ADCC = np.array([mean_3a, mean_3b])

    with M:
        x_x = pm.Lognormal("activity_scores", sigma=1.0, shape=(24, 2))

        pm.Normal("scale", mu=1.0, sigma=0.1, observed=T.sum(x_x, axis=1)) # Enforce unit scaled loadings

        residuals = mean_ADCC - T.dot(A_antiD, x_x).T
        sd = T.minimum(T.std(residuals), 1.0)  # Add bounds for the stderr to help force the fitting solution
        pm.Normal("fit", sigma=sd, observed=residuals)

    trace = pm.sample(2000, init="advi+adapt_diag", model=M, return_inferencedata=True, target_accept=0.9, chains=1)

    az.plot_forest(trace)