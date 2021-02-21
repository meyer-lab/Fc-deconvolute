from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x

def makeFigure():
    """Check the fit using infered x values"""
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, glycan_list = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    glycans_3a = infer_x(A_antiD, mean_3a)
    glycans_3b = infer_x(A_antiD, mean_3b)

    infer_adcc_3a = A_antiD @ glycans_3a
    infer_adcc_3b = glycans_3b.dot(A_antiD)

    ax[0].scatter(glycan_list, infer_adcc_3a)
    ax[0].set_title("ADCC (Fig. 3A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)

    ax[1].scatter(glycan_list, infer_adcc_3b)
    ax[1].set_title("ADCC (Fig. 3B)")
    ax[1].set_xlabel("Glycans")
    ax[1].set_xticklabels(ax[1].get_xticks(), rotation=90)

    # Add subplot labels
    subplotLabel(ax)

    return f
