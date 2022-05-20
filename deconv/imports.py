import pandas as pd
from sklearn.preprocessing import scale


def load_dekkers(mean=True):
    results = dict()

    # anti-D
    antiD = pd.read_csv("./deconv/data/anti-D.csv")
    results["antiD"] = antiD.iloc[:, 7:].values / 100.0
    results["glycans"] = list(antiD.columns.values[7:])
    results["mixtures"] = antiD.iloc[:, 0]

    # load profiling data
    profiling = pd.read_csv("./deconv/data/dekkers.csv")

    if mean:
        profiling = profiling.groupby(["index", "receptor"]).mean().reset_index()
        profiling = profiling.pivot(index="index", columns="receptor", values="binding")
        scale(profiling, with_mean=False, copy=False, axis=1)

    results["profiling"] = profiling
    return results
