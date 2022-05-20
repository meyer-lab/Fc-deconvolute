import pandas as pd


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
        profiling_std = profiling.groupby(["index", "receptor"]).sem().reset_index()
        profiling = profiling.groupby(["index", "receptor"]).mean().reset_index()
        profiling = profiling.pivot(index="index", columns="receptor", values="binding")
        profiling_std = profiling_std.pivot(index="index", columns="receptor", values="binding")

        col_means = profiling.mean()
        results["profiling_std"] = profiling_std / col_means
        profiling = profiling / col_means

    print(results["profiling_std"].mean())

    results["profiling"] = profiling
    return results
