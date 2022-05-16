
from .imports import load_dekkers
import numpy as np
import pymc as pm
import aesara.tensor as T


def getEmceeTrace():
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]
    A_antiD /= np.sum(A_antiD, axis=1, keepdims=True)

    df = data_dekkers["profiling"]
    data = df.groupby(["index", "receptor"]).mean().reset_index()
    res = data.pivot(index="receptor", columns="index", values="binding").values
    res /= np.mean(res, axis=1, keepdims=True)
    res = res.T

    M = pm.Model()

    with M:
        activity = pm.Lognormal("activity", sigma=1.0, mu=-1.0, shape=(24, 12))
        predict = T.dot(A_antiD, activity)
        pm.Lognormal("fit", mu=predict, sigma=0.1, observed=res)

    trace = pm.sample(3000, tune=3000, model=M, return_inferencedata=True, target_accept=0.95)
    return trace
