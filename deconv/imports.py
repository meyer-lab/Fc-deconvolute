import pandas as pd
from scipy.optimize import nnls


def load_tables():
    antiD = pd.read_csv("./deconv/data/anti-D.csv")
    antiT = pd.read_csv("./deconv/data/anti-TNP.csv")
    glycan_list = list(antiD.columns.values[7:])
    A_antiD = antiD.iloc[:, 7:].values
    A_antiTNP = antiT.iloc[:, 7:].values
    return (A_antiD, A_antiTNP, glycan_list)


def load_figures():
    # binding in doners with specific r3a gamma receptor
    figA = pd.read_csv("./deconv/data/fig3a.csv")
    # combination of 4 independent receptors?
    # number of the mixture (4 data points for each mixture), cytotoxicity percentage
    figB = pd.read_csv("./deconv/data/fig3b.csv")
    adcc_3a = figA.iloc[:, 0]
    adcc_3b = figB.iloc[:, 0]
    return (adcc_3a, adcc_3b)


def infer_x(A, adcc):
    return nnls(A, adcc, maxiter=None)[0]
