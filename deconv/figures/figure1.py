import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x, infer_x_fixed, ADCC_groups
from sklearn.utils import resample


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, glycan_list, _ = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    setGroup = ADCC_groups()

    glycans_3a = infer_x_fixed(A_antiD, mean_3a, setGroup)
    glycans_3b = infer_x_fixed(A_antiD, mean_3b, setGroup)

    ax[0].bar(glycan_list, glycans_3a)
    ax[0].set_title("ADCC (Fig. 3A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(glycan_list, rotation=90)

    ax[1].bar(glycan_list, glycans_3b)
    ax[1].set_title("ADCC (Fig. 3B)")
    ax[1].set_xlabel("Glycans")
    ax[1].set_xticklabels(glycan_list, rotation=90)

    # Add subplot labels
    subplotLabel(ax)

    return f

def makeFigureWithError():
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, glycan_list = load_tables()
    adcc_3a, adcc_3b = load_figures()

    num_iters = 100

    glycans_3a = []
    glycans_3b = []

    for i in range(num_iters):
        new3a = resample(adcc_3a, n_samples=40, replace=True, stratify=adcc_3a.index)
        mean_3a = new3a.groupby(level=0).sum() / new3a.groupby(level=0).count()

        new3b = resample(adcc_3b, n_samples=40, replace=True, stratify=adcc_3b.index)
        mean_3b = new3b.groupby(level=0).sum() / new3b.groupby(level=0).count()

        glycans_3a.append(infer_x(A_antiD, mean_3a))
        glycans_3b.append(infer_x(A_antiD, mean_3b))
    
    glycans3a = np.array(glycans_3a)
    mean3a = np.mean(glycans3a, axis=0)
    error3a = np.std(glycans3a, axis=0)

    glycans3b = np.array(glycans_3b)
    mean3b = np.mean(glycans3b, axis=0)
    error3b = np.std(glycans3b, axis=0)

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