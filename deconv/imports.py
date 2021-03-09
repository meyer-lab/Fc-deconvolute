import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.optimize import least_squares


def load_tables():
    antiD = pd.read_csv("./deconv/data/anti-D.csv")
    antiT = pd.read_csv("./deconv/data/anti-TNP.csv")
    glycan_list = list(antiD.columns.values[7:])
    mixtures = antiD.iloc[:,0]
    A_antiD = antiD.iloc[:, 7:].values
    A_antiTNP = antiT.iloc[:, 7:].values
    return (A_antiD, A_antiTNP, glycan_list, mixtures)


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


def cost(pIn, X, y, setGroups, norm=False):
    outt = X @ pIn[setGroups] - y
    if norm:
        outt = np.linalg.norm(outt)
    return outt


def infer_x_fixed(X, y, setGroups, numGroups=None, retP=False):
    assert setGroups.dtype == np.int
    assert setGroups.ndim == 1
    assert y.ndim == 1
    assert y.size == X.shape[0]
    assert X.shape[1] == setGroups.size

    if numGroups is None:
        numGroups = len(np.unique(setGroups))

    res = least_squares(lambda pp: cost(pp, X, y, setGroups), np.ones(numGroups), ftol=1e-9, bounds=(0, np.inf))
    assert res.success

    if retP:
        return res.x

    return res.x[setGroups]


def infer_x_EM(X, y, nGroups):
    """ Sets up a strategy to fit the levels if we don't know them. """
    setGroups = np.random.choice(nGroups, size=X.shape[1])
    assert X.shape[1] == setGroups.size
    assert y.size == X.shape[0]
    pIn = infer_x_fixed(X, y, setGroups, numGroups=nGroups, retP=True)
    baseCost = cost(pIn, X, y, setGroups, norm=True)

    for _ in range(1000):
        pos = np.random.choice(setGroups.size, size=1, replace=False)
        new = np.random.choice(np.max(setGroups), size=1)
        newGroups = np.copy(setGroups)
        newGroups[pos] = new

        pNew = infer_x_fixed(X, y, newGroups, numGroups=nGroups, retP=True)

        newCost = cost(pNew, X, y, newGroups, norm=True)
        if newCost < baseCost:
            setGroups = newGroups
            pIn = pNew
            baseCost = newCost

    return pIn[setGroups]
