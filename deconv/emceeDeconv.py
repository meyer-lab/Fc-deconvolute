import numpy as np
from theano import tensor as T
from .imports import load_dekkers
import pymc3 as pm


def getEmceeTrace():

    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]

    mean_binding = list(data_dekkers["bindings"])
    mean_binding = np.array([m.groupby(level=0).mean() for m in mean_binding])

    mean_3a = data_dekkers["meanADCC3a"]
    mean_3b = data_dekkers["meanADCC3b"]

    mean_ADCC = np.array([mean_3a, mean_3b])

    mean_4a = data_dekkers["meanCompAct4a"]
    mean_4b = data_dekkers["meanCompAct4b"]
    mean_Comp = np.array([mean_4a, mean_4b])

    res = np.concatenate((mean_binding, mean_ADCC, mean_Comp))

    M = pm.Model()

    with M:
        x_x = pm.Lognormal("activity_scores", sigma=1.0, shape=(24, 3))
        x_y = pm.Lognormal("activity_loadings", sigma=1.0, shape=(3, 12))

        pm.Normal("scale", mu=1.0, sigma=0.1, observed=T.sum(x_y, axis=1))  # Enforce unit scaled loadings

        residuals = res - T.dot(T.dot(A_antiD, x_x), x_y).T
        sd = T.minimum(T.std(residuals), 1.0)  # Add bounds for the stderr to help force the fitting solution
        pm.Normal("fit", sigma=sd, observed=residuals)

    trace = pm.sample(2000, init="advi+adapt_diag", model=M, return_inferencedata=True, target_accept=0.9, chains=1)

    return trace
