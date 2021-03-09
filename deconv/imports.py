import pandas as pd
from scipy.optimize import nnls


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


def load_bindingData():
    fig_2A = pd.read_csv("./deconv/data/Fig2A-FcgRI.csv")
    fig_2B = pd.read_csv("./deconv/data/Fig2B-FcgRIIa-131H.csv")
    fig_2C = pd.read_csv("./deconv/data/Fig2C-FcgRIIa-131R.csv")
    fig_2D = pd.read_csv("./deconv/data/Fig2D-FcgRIIb.csv")
    fig_2E = pd.read_csv("./deconv/data/Fig2E-FcgRIIIa-158F.csv")
    fig_2F = pd.read_csv("./deconv/data/Fig2F-FcgRIIIa-158V.csv")
    fig_2G = pd.read_csv("./deconv/data/Fig2G-FcgRIIIb-NA1.csv")
    fig_2H = pd.read_csv("./deconv/data/Fig2H-FcgRIIIb-NA2.csv")
    fig2A = fig_2A.iloc[:, :].values
    fig2B = fig_2B.iloc[:, :].values
    fig2C = fig_2C.iloc[:, :].values
    fig2D = fig_2D.iloc[:, :].values
    fig2E = fig_2E.iloc[:, :].values
    fig2F = fig_2F.iloc[:, :].values
    fig2G = fig_2G.iloc[:, :].values
    fig2H = fig_2H.iloc[:, :].values
    return(fig2A, fig2B, fig2C, fig2D, fig2E, fig2F, fig2G, fig2H)
