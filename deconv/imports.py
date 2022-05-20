import numpy as np
import pandas as pd


def load_dekkers(mean=True):
    results = dict()

    # anti-D
    antiD = pd.read_csv("./deconv/data/anti-D.csv")
    results["antiD"] = antiD.iloc[:, 7:].values
    results["glycans"] = list(antiD.columns.values[7:])
    results["mixtures"] = antiD.iloc[:, 0]

    # load profiling data
    profiling = pd.read_csv("./deconv/data/dekkers.csv")

    if mean:
        profilingG = profiling.groupby(["index", "receptor"]).mean()
        profilingG["binding"] = profilingG.transform(lambda x: x / np.std(x))  # Dividing by std per experiement
        profiling = profilingG.reset_index()

    results["profiling"] = profiling
    return results
