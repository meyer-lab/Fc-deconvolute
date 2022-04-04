from theano import tensor as T
from .imports import load_dekkers
import pymc3 as pm


def getEmceeTrace():
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]

    df = data_dekkers["profiling"]
    data = df.groupby(["index", "receptor"]).mean().reset_index()
    res = data.pivot(index="receptor", columns="index", values="binding").values

    M = pm.Model()

    with M:
        x_x = pm.Lognormal("activity_scores", sigma=1.0, shape=(24, 3))
        x_y = pm.Lognormal("activity_loadings", sigma=1.0, shape=(3, 12))

        pm.Normal("scale", mu=1.0, sigma=0.1, observed=T.sum(x_y, axis=1))  # Enforce unit scaled loadings

        residuals = res - T.dot(T.dot(A_antiD, x_x), x_y).T
        sd = T.minimum(T.std(residuals), 1.0)  # Add bounds for the stderr to help force the fitting solution
        pm.Normal("fit", sigma=sd, observed=residuals)

    trace = pm.sample(2000, init="advi+adapt_diag", model=M, return_inferencedata=True, target_accept=0.95, chains=None, cores=4)

    return trace
