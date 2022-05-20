
from .imports import load_dekkers
import pymc as pm
from aesara.tensor import dot


def getEmceeTrace():
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]
    res = data_dekkers["profiling"].values

    M = pm.Model()

    with M:
        activity = pm.Lognormal("activity", sigma=1.0, shape=(24, 12))
        predict = dot(A_antiD, activity).T
        pm.Lognormal("fit", mu=predict, sigma=0.2, observed=res)

    trace = pm.sample(500, tune=5000, model=M, return_inferencedata=True, target_accept=0.95)
    return trace
