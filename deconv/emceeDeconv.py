from .imports import load_dekkers
import pymc as pm
import aesara.tensor as T


def getEmceeTrace():
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]
    res = data_dekkers["profiling"].values

    M = pm.Model()

    with M:
        activity = pm.TruncatedNormal("activity", mu=0.0, sigma=5.0, lower=0.0, shape=(24, 12))
        predict = T.dot(A_antiD, activity)
        pm.Normal("fit", mu=predict, sigma=0.1, observed=res)

    return pm.sample(model=M, return_inferencedata=True, target_accept=0.95)
