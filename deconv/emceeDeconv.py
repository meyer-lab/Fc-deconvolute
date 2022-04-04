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
        activity = pm.Lognormal("activity", sigma=1.0, shape=(24, 12))
        predict = T.dot(A_antiD, activity).T
        pm.Lognormal("fit", mu=predict, sigma=0.5, observed=res)

    trace = pm.sample(2000, model=M, return_inferencedata=True, target_accept=0.95)
    return trace
