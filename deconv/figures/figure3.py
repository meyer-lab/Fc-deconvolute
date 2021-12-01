import numpy as np
from .common import subplotLabel, getSetup
from ..imports import infer_x, load_dekkers


def makeFigure():
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    data_dekkers = load_dekkers()

    A_antiD, A_antiTNP, glycan_list = data_dekkers["antiD"], data_dekkers["antiTNP"], data_dekkers["glycans"]

    mean_3a = data_dekkers["meanADCC3a"]
    mean_3b = data_dekkers["meanADCC3b"]

    A = np.concatenate((A_antiD, A_antiTNP), axis=0)

    double_3a = np.concatenate((mean_3a, mean_3a), axis=0)
    double_3b = np.concatenate((mean_3b, mean_3b), axis=0)

    glycans_3a = infer_x(A, double_3a)
    glycans_3b = infer_x(A, double_3b)

    ind = np.arange(len(glycan_list))

    ax[0].bar(ind, glycans_3a, label='Anti-D and Anti-TNP')
    ax[0].set_title("Anti-D and Anti-TNP (Fig. 3A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(glycan_list, rotation=90)
    ax[0].set_xticks(ind)
    ax[0].legend()

    ax[1].bar(ind, glycans_3b, label='Anti-D and Anti-TNP')
    ax[1].set_title("Anti-D and Anti-TNP (Fig. 3B)")
    ax[1].set_xlabel("Glycans")
    ax[1].set_xticklabels(glycan_list, rotation=90)
    ax[1].set_xticks(ind)
    ax[1].legend()

    # Add subplot labels
    subplotLabel(ax)

    return f
