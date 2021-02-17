import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x

def makeFigure():
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, A_antiTNP, glycan_list = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    A = np.concatenate((A_antiD, A_antiTNP), axis=0)

    double_3a = np.concatenate((mean_3a, mean_3a), axis=0)
    double_3b = np.concatenate((mean_3b, mean_3b), axis=0)
    double_antiD = np.concatenate((A_antiD, A_antiD), axis=0)

    glycans_3a = infer_x(A, double_3a)
    glycans_3b = infer_x(A, double_3b)

    original_3a = infer_x(double_antiD, double_3a)
    original_3b = infer_x(double_antiD, double_3b)

    width = 0.35
    ind = np.arange(len(glycan_list))

    ax[0].bar(ind - width/2, glycans_3a, width, label = 'Anti-D and Anti-TNP')
    ax[0].bar(ind + width/2, original_3a, width, label = 'Double Anti-D')
    ax[0].set_title("Anti-D and Anti-TNP (Fig. 3A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(glycan_list, Rotation=90)
    ax[0].set_xticks(ind)
    ax[0].legend()

    ax[1].bar(ind - width/2, glycans_3b, width, label = 'Anti-D and Anti-TNP')
    ax[1].bar(ind + width/2, original_3b, width, label = 'Double Anti-D')
    ax[1].set_title("Anti-D and Anti-TNP (Fig. 3B)")
    ax[1].set_xlabel("Glycans")
    ax[1].set_xticklabels(glycan_list, Rotation=90)
    ax[1].set_xticks(ind)
    ax[1].legend()

    # Add subplot labels
    subplotLabel(ax)

    return f
