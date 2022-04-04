import pandas as pd


def load_dekkers():
    results = dict()

    # anti-D
    antiD = pd.read_csv("./deconv/data/anti-D.csv")
    results["antiD"] = antiD.iloc[:, 7:].values
    results["glycans"] = list(antiD.columns.values[7:])
    results["mixtures"] = antiD.iloc[:, 0]

    # load profiling data
    profiling = pd.read_csv("./deconv/data/dekkers.csv")
    results["profiling"] = profiling
    return results
