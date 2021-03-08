import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, infer_x, load_bindingData

def makeFigure():
    ax, f = getSetup((6, 12), (4, 2))

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
    glycans_2C = infer_x(A_antiD, mean_2C)
    glycans_2D = infer_x(A_antiD, mean_2D)
    glycans_2E = infer_x(A_antiD, mean_2E)
    glycans_2F = infer_x(A_antiD, mean_2F)
    glycans_2G = infer_x(A_antiD, mean_2G)
    glycans_2H = infer_x(A_antiD, mean_2H)

    ax[0].bar(glycan_list, glycans_2A)
    ax[0].set_title("FcgRI Binding (Fig. 2A)")

    ax[1].bar(glycan_list, glycans_2B)
    ax[1].set_title("FcgRIIa-131H Binding (Fig. 2B)")

    ax[2].bar(glycan_list, glycans_2C)
    ax[2].set_title("FcgRIIa-131R Binding (Fig. 2C)")

    ax[3].bar(glycan_list, glycans_2D)
    ax[3].set_title("FcgRIIb Binding (Fig. 2D)")

    ax[4].bar(glycan_list, glycans_2E)
    ax[4].set_title("FcgRIIIa-158F Binding (Fig. 2E)")

    ax[5].bar(glycan_list, glycans_2F)
    ax[5].set_title("FcgRIIIa-158V Binding (Fig. 2F)")

    ax[6].bar(glycan_list, glycans_2G)
    ax[6].set_title("FcgRIIIb-NA1 Binding (Fig. 2G)")

    ax[7].bar(glycan_list, glycans_2H)
    ax[7].set_title("FcgRIIIb-NA2 Binding (Fig. 2H)")

    for x in range(8):
        ax[x].set_xlabel("Glycans")
        ax[x].set_xticklabels(glycan_list, rotation=90)

    # Add subplot labels
    subplotLabel(ax)

    return f
