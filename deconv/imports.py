import numpy as np
import pandas as pd


def load_dekkers():
    results = dict()

    # anti d anti tnp
    antiD = pd.read_csv("./deconv/data/anti-D.csv")
    antiT = pd.read_csv("./deconv/data/anti-TNP.csv")
    glycan_list = list(antiD.columns.values[7:])
    mixtures = antiD.iloc[:, 0]
    A_antiD = antiD.iloc[:, 7:].values
    A_antiTNP = antiT.iloc[:, 7:].values
    results["antiD"] = A_antiD
    results["antiTNP"] = A_antiTNP
    results["glycans"] = glycan_list
    results["mixtures"] = mixtures

    # load figure 3
    figA = pd.read_csv("./deconv/data/fig3a.csv")
    figB = pd.read_csv("./deconv/data/fig3b.csv")
    adcc_3a = figA.iloc[:, 0]
    adcc_3b = figB.iloc[:, 0]
    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4
    results["meanADCC3a"] = mean_3a
    results["meanADCC3b"] = mean_3b

    # load figure 4
    figA = pd.read_csv('./deconv/data/fig4a c1q.csv')
    figB = pd.read_csv('./deconv/data/fig4b c4.csv')
    complement_4a = figA.to_numpy()[:, 1:]
    complement_4b = figB.to_numpy()[:, 1:]
    mean_4a = np.mean(complement_4a, axis=1)
    mean_4b = np.mean(complement_4b, axis=1)
    results["meanCompAct4a"] = mean_4a
    results["meanCompAct4b"] = mean_4b

    # load figure 2
    fig2A = pd.read_csv("./deconv/data/Fig2A-FcgRI.csv", index_col=0)
    RI = fig2A.iloc[:, 1]
    fig2B = pd.read_csv("./deconv/data/Fig2B-FcgRIIa-131H.csv", index_col=0)
    RIIa_131H = fig2B.iloc[:, 1]
    fig2C = pd.read_csv("./deconv/data/Fig2C-FcgRIIa-131R.csv", index_col=0)
    RIIa_131R = fig2C.iloc[:, 1]
    fig2D = pd.read_csv("./deconv/data/Fig2D-FcgRIIbc.csv", index_col=0)
    RIIb = fig2D.iloc[:, 1]
    fig2E = pd.read_csv("./deconv/data/Fig2E-FcgRIIIa-158F.csv", index_col=0)
    RIIIa_158F = fig2E.iloc[:, 1]
    fig2F = pd.read_csv("./deconv/data/Fig2F-FcgRIIIa-158V.csv", index_col=0)
    RIIIa_158V = fig2F.iloc[:, 1]
    fig2G = pd.read_csv("./deconv/data/Fig2G-FcgRIIIb-NA1.csv", index_col=0)
    RIIIb_NA1 = fig2G.iloc[:, 1]
    fig2H = pd.read_csv("./deconv/data/Fig2H-FcgRIIIb-NA2.csv", index_col=0)
    RIIIb_NA2 = fig2H.iloc[:, 1]
    binding = [RI, RIIa_131H, RIIa_131R, RIIb, RIIIa_158F, RIIIa_158V, RIIIb_NA1, RIIIb_NA2]
    results["bindings"] = binding

    return (results)
