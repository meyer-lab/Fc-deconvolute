import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x, infer_x_fixed


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, glycan_list, _ = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    setGroup = np.zeros(A_antiD.shape[1], dtype=np.int)
    setGroup[[13, 18]] = 1
    setGroup[[14, 19]] = 2
    setGroup[[15, 20]] = 3
    setGroup[[16, 21]] = 4
    setGroup[[17, 22]] = 4
    setGroup[23] = 5
    setGroup[12] = 6

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
