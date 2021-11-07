import numpy as np
import pandas as pd
from .common import subplotLabel, getSetup
from ..imports import infer_x_fixed, ADCC_groups, load_dekkers
from sklearn.utils import resample


def makeFigure():
    """ Look at bootstrapping of ADCC. """
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    data_dekkers = load_dekkers()
    
    A_antiD, glycan_list = data_dekkers["antiD"],  data_dekkers["glycans"]

    mean_3a = data_dekkers["meanADCC3a"]
    mean_3b = data_dekkers["meanADCC3b"]

    figA = pd.read_csv("./deconv/data/fig3a.csv")
    figB = pd.read_csv("./deconv/data/fig3b.csv")
    adcc_3a = figA.iloc[:, 0]
    adcc_3b = figB.iloc[:, 0]

    setGroup = ADCC_groups()

    num_iters = 100
    mean3a = infer_x_fixed(A_antiD, mean_3a, setGroup)
    mean3b = infer_x_fixed(A_antiD, mean_3b, setGroup)

    glycans_3a = []
    glycans_3b = []

    for i in range(num_iters):
        new3a = resample(adcc_3a, n_samples=40, replace=True, stratify=adcc_3a.index)
        mean_3a = new3a.groupby(level=0).sum() / new3a.groupby(level=0).count()

        new3b = resample(adcc_3b, n_samples=40, replace=True, stratify=adcc_3b.index)
        mean_3b = new3b.groupby(level=0).sum() / new3b.groupby(level=0).count()

        glycans_3a.append(infer_x_fixed(A_antiD, mean_3a, setGroup))
        glycans_3b.append(infer_x_fixed(A_antiD, mean_3b, setGroup))

    glycans3a = np.array(glycans_3a)
    error3a = (mean3a - np.quantile(glycans3a, 0.33, axis=0), np.quantile(glycans3a, 0.67, axis=0) - mean3a)

    glycans3b = np.array(glycans_3b)
    error3b = (mean3b - np.quantile(glycans3b, 0.33, axis=0), np.quantile(glycans3b, 0.67, axis=0) - mean3b)

    ax[0].bar(glycan_list, mean3a, yerr=error3a)
    ax[0].set_title("ADCC (Fig. 3A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(glycan_list, rotation=90)

    ax[1].bar(glycan_list, mean3b, yerr=error3b)
    ax[1].set_title("ADCC (Fig. 3B)")
    ax[1].set_xlabel("Glycans")
    ax[1].set_xticklabels(glycan_list, rotation=90)

    # Add subplot labels
    subplotLabel(ax)

    return f
