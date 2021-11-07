from .common import subplotLabel, getSetup
from ..imports import infer_x_fixed, ADCC_groups, load_dekkers


def makeFigure():
    """Get a list of the axis objects and create a figure"""
    # Get list of axis objects
    ax, f = getSetup((6, 3), (1, 2))

    data_dekkers = load_dekkers()

    A_antiD, glycan_list = data_dekkers["antiD"], data_dekkers["glycans"]

    mean_3a = data_dekkers["meanADCC3a"]
    mean_3b = data_dekkers["meanADCC3b"]

    setGroup = ADCC_groups()

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
