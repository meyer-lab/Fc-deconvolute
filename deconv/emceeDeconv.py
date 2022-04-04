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
        activity = pm.Lognormal("activity", sigma=1.0, shape=(24, 12))
        predict = T.dot(A_antiD, activity).T
        pm.Lognormal("fit", mu=predict, sigma=0.2, observed=res)

    trace = pm.sample(200, model=M, return_inferencedata=True, target_accept=0.95)
    return trace
