import numpy as np
from .common import subplotLabel, getSetup
from ..imports import infer_x_fixed, ADCC_groups, load_dekkers


def makeFigure():
    """Check the fit using infered x values"""
    ax, f = getSetup((6, 3), (1, 2))

    data_dekkers = load_dekkers()

    A_antiD, A_antiTNP, mixtures = data_dekkers["antiD"], data_dekkers["antiTNP"], data_dekkers["mixtures"]

    mean_3a = data_dekkers["meanADCC3a"]
    mean_3b = data_dekkers["meanADCC3b"]

    double_3a = np.concatenate((mean_3a, mean_3a), axis=0)
    double_3b = np.concatenate((mean_3b, mean_3b), axis=0)

    A = np.concatenate((A_antiD, A_antiTNP), axis=0)

    setGroup = ADCC_groups()

    glycans_3a = infer_x_fixed(A, double_3a, setGroup)
    glycans_3b = infer_x_fixed(A, double_3b, setGroup)

    infer_adcc_3a = A_antiD @ glycans_3a
    infer_adcc_3b = A_antiD @ glycans_3b

    width = 0.35
    ind = np.arange(len(mixtures))

    ax[0].bar(ind - width / 2, infer_adcc_3a, width, label='Inferred ADCC')
    ax[0].bar(ind + width / 2, mean_3a, width, label='Original ADCC')
    ax[0].set_title("Original and Inferred ADCC (Fig. 3A)")
    ax[0].set_xlabel("Mixtures")
    ax[0].set_xticklabels(mixtures, rotation=90)
    ax[0].set_xticks(ind)
    ax[0].legend()

    ax[1].bar(ind - width / 2, infer_adcc_3b, width, label='Inferred ADCC')
    ax[1].bar(ind + width / 2, mean_3b, width, label='Original ADCC')
    ax[1].set_title("Original and Inferred ADCC (Fig. 3B)")
    ax[1].set_xlabel("Mixtures")
    ax[1].set_xticklabels(mixtures, rotation=90)
    ax[1].set_xticks(ind)
    ax[1].legend()

    # Add subplot labels
    subplotLabel(ax)

    return f
