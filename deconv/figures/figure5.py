import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, infer_x, load_bindingData

def makeFigure():
    ax, f = getSetup((6, 3), (2, 4))

    fig2A, fig2B, fig2C, fig2D, fig2E, fig2F, fig2G, fig2H = load_bindingData()

    A_antiD, _, glycan_list, mixtures = load_tables()

    mean_2A = []
    mean_2B = []
    mean_2C = []
    mean_2D = []
    mean_2E = []
    mean_2F = []
    mean_2G = []
    mean_2H = []

    for x in np.arange(len(mixtures)):
        mean_2A.append(np.nanmean(fig2A[:, x]))
        mean_2B.append(np.nanmean(fig2B[:, x]))
        mean_2C.append(np.nanmean(fig2C[:, x]))
        mean_2D.append(np.nanmean(fig2D[:, x]))
        mean_2E.append(np.nanmean(fig2E[:, x]))
        mean_2F.append(np.nanmean(fig2F[:, x]))
        mean_2G.append(np.nanmean(fig2G[:, x]))
        mean_2H.append(np.nanmean(fig2H[:, x]))

    glycans_2A = infer_x(A_antiD, mean_2A)
    glycans_2B = infer_x(A_antiD, mean_2B)

    ax[0].bar(glycan_list, glycans_2A)
    ax[0].set_title("FcgRI Binding (Fig. 2A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(glycan_list, rotation=90)

    ax[1].bar(glycan_list, glycans_2B)
    ax[1].set_title("FcgRIIa Binding (Fig. 2B)")
    ax[1].set_xlabel("Glycans")
    ax[1].set_xticklabels(glycan_list, rotation=90)

    # Add subplot labels
    subplotLabel(ax)

    return f
