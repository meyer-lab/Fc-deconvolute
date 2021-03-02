import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x

def makeFigure():
    """Artificially alter G1S"""
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, _, mixtures = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    glycans_3a = infer_x(A_antiD, mean_3a)
    glycans_3b = infer_x(A_antiD, mean_3b)

    original_3a = A_antiD @ glycans_3a
    original_3b = A_antiD @ glycans_3b

    glycans_3a[18] = glycans_3a[12:23].sum() / 12
    glycans_3b[18] = glycans_3b[12:23].sum() / 12

    g1s_3a = A_antiD @ glycans_3a
    g1s_3b = A_antiD @ glycans_3b

    width = 0.35
    ind = np.arange(len(mixtures))

    ax[0].bar(ind - width/2, g1s_3a, width, label = 'Altered G1S')
    ax[0].bar(ind + width/2, original_3a, width, label = 'Original G1S')
    ax[0].set_title("Original and Altered G1S (Fig. 3A)")
    ax[0].set_xlabel("Mixtures")
    ax[0].set_xticklabels(mixtures, Rotation=90)
    ax[0].set_xticks(ind)
    ax[0].legend()

    ax[1].bar(ind - width/2, g1s_3b, width, label = 'Altered G1S')
    ax[1].bar(ind + width/2, original_3b, width, label = 'Original G1S')
    ax[1].set_title("Original and Altered G1S (Fig. 3B)")
    ax[1].set_xlabel("Mixtures")
    ax[1].set_xticklabels(mixtures, Rotation=90)
    ax[1].set_xticks(ind)
    ax[1].legend()

    # Add subplot labels
    subplotLabel(ax)

    return f