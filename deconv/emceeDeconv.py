
from .imports import load_dekkers
import pymc as pm
from aesara.tensor import dot


def getEmceeTrace():
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]

    df = data_dekkers["profiling"]
    effector_responses = df["receptor"].unique()

    # normalization
    for res in effector_responses:
        for i in range(1,21):
            mask = (df["receptor"] == res) & (df["index"] == i)
            df1 = df.loc[mask]
            df1["binding"] /= df1["binding"].std()
            df.loc[mask] = df1

    data = df.groupby(["index", "receptor"]).mean().reset_index()

    res = data.pivot(index="receptor", columns="index", values="binding").values

    M = pm.Model()

    with M:
        activity = pm.Lognormal("activity", sigma=1.0, shape=(24, 12))
        predict = dot(A_antiD, activity).T
        pm.Lognormal("fit", mu=predict, sigma=0.2, observed=res)

    trace = pm.sample(500, tune=5000, model=M, return_inferencedata=True, target_accept=0.95)
    return trace
