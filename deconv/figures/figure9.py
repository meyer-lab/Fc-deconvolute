import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, load_bindingData, infer_x, infer_x_fixed, ADCC_groups, R1_groups

def makeFigure():
    """Binding fit"""
    ax, f = getSetup((6, 12), (4, 2))
    A_antiD, _, _, mixtures = load_tables()
    mean_binding = list(load_bindingData())
    R3_groups = ADCC_groups()
    R1_group = R1_groups()

    for ii, data in enumerate(mean_binding):
        mean_binding[ii] = mean_binding[ii].groupby(level=0).sum() / mean_binding[ii].groupby(level=0).count()

    deconvA_grouped = infer_x_fixed(A_antiD, mean_binding[0], R1_group)
    deconvB = infer_x(A_antiD, mean_binding[1])
    deconvC = infer_x(A_antiD, mean_binding[2])
    deconvD = infer_x(A_antiD, mean_binding[3])
    deconvE_grouped = infer_x_fixed(A_antiD, mean_binding[4], R3_groups)
    deconvF_grouped = infer_x_fixed(A_antiD, mean_binding[5], R3_groups)
    deconvG_grouped = infer_x_fixed(A_antiD, mean_binding[6], R3_groups)
    deconvH_grouped = infer_x_fixed(A_antiD, mean_binding[7], R3_groups)
    deconv = [deconvA_grouped, deconvB, deconvC, deconvD, deconvE_grouped, deconvF_grouped, deconvG_grouped, deconvH_grouped]

    width = 0.35
    ind = np.arange(len(mixtures))

    for x in range(8):
        infer_binding = A_antiD @ deconv[x]
        ax[x].bar(ind - width/2, infer_binding, width, label = 'Inferred')
        ax[x].bar(ind + width/2, mean_binding[x], width, label = 'Measured')
        ax[x].set_xlabel("Mixtures")
        ax[x].set_xticklabels(mixtures, rotation=90)
        ax[x].set_ylabel("Receptor Binding")
        ax[x].set_xticks(ind)
        ax[x].legend()

    ax[0].set_title("FcγRI")
    ax[1].set_title("FcγRIIa-131H")
    ax[2].set_title("FcγRIIa-131R")
    ax[3].set_title("FcγRIIb/c")
    ax[4].set_title("FcγRIIIa-158F")
    ax[5].set_title("FcγRIIIa-158V")
    ax[6].set_title("FcγRIIIb-NA1")
    ax[7].set_title("FcγRIIIb-NA2")

    subplotLabel(ax)

    return f
