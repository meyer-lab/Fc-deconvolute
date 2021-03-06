import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x

def makeFigure():
    """Check the fit using infered x values"""
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, _, mixtures = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    glycans_3a = infer_x(A_antiD, mean_3a)
    glycans_3b = infer_x(A_antiD, mean_3b)

    infer_adcc_3a = A_antiD @ glycans_3a
    infer_adcc_3b = A_antiD @ glycans_3b

    width = 0.35
    ind = np.arange(len(mixtures))

    ax[0].bar(ind - width/2, infer_adcc_3a, width, label = 'Inferred ADCC')
    ax[0].bar(ind + width/2, mean_3a, width, label = 'Original ADCC')
    ax[0].set_title("Original and Inferred ADCC (Fig. 3A)")
    ax[0].set_xlabel("Mixtures")
    ax[0].set_xticklabels(mixtures, Rotation=90)
    ax[0].set_xticks(ind)
    ax[0].legend()

    ax[1].bar(ind - width/2, infer_adcc_3b, width, label = 'Inferred ADCC')
    ax[1].bar(ind + width/2, mean_3b, width, label = 'Original ADCC')
    ax[1].set_title("Original and Inferred ADCC (Fig. 3B)")
    ax[1].set_xlabel("Mixtures")
    ax[1].set_xticklabels(mixtures, Rotation=90)
    ax[1].set_xticks(ind)
    ax[1].legend()

    # Add subplot labels
    subplotLabel(ax)

    return f
