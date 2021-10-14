import numpy as np
from .common import subplotLabel, getSetup
from ..imports import load_tables, infer_x, load_bindingData


def makeFigure():
    ax, f = getSetup((6, 12), (4, 2))

    binding_data = list(load_bindingData())

    A_antiD, _, glycan_list, _ = load_tables()

    for ii, data in enumerate(binding_data):
        binding_data[ii] = np.nanmean(data, axis=0)

    deconv = [infer_x(A_antiD, dd) for dd in binding_data]

    ax[0].bar(glycan_list, deconv[0])
    ax[0].set_title("FcgRI Binding (Fig. 2A)")

    ax[1].bar(glycan_list, deconv[1])
    ax[1].set_title("FcgRIIa-131H Binding (Fig. 2B)")

    ax[2].bar(glycan_list, deconv[2])
    ax[2].set_title("FcgRIIa-131R Binding (Fig. 2C)")

    ax[3].bar(glycan_list, deconv[3])
    ax[3].set_title("FcgRIIb Binding (Fig. 2D)")

    ax[4].bar(glycan_list, deconv[4])
    ax[4].set_title("FcgRIIIa-158F Binding (Fig. 2E)")

    ax[5].bar(glycan_list, deconv[5])
    ax[5].set_title("FcgRIIIa-158V Binding (Fig. 2F)")

    ax[6].bar(glycan_list, deconv[6])
    ax[6].set_title("FcgRIIIb-NA1 Binding (Fig. 2G)")

    ax[7].bar(glycan_list, deconv[7])
    ax[7].set_title("FcgRIIIb-NA2 Binding (Fig. 2H)")

    for x in range(8):
        ax[x].set_xlabel("Glycans")
        ax[x].set_xticklabels(glycan_list, rotation=90)

    # Add subplot labels
    subplotLabel(ax)

    return f
