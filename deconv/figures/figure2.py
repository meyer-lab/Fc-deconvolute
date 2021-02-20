from .common import subplotLabel, getSetup
from ..imports import load_tables, load_figures, infer_x
import numpy as np 
import pandas as pd

def fitAssess():
    """Check the fit using infered x values"""
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, _ = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    glycans_3a = infer_x(A_antiD, mean_3a)
    glycans_3b = infer_x(A_antiD, mean_3b)

    infer_adcc_3a = np.matmul(A_antiD, glycans_3a)
    infer_adcc_3b = np.matmul(A_antiD, glycans_3b)

    mixture = list(range(1,21))
    width = 0.35
    ind = np.arange(len(mixture))

    ax[0].bar(ind - width/2, infer_adcc_3a, width, label = 'Inferred ADCC Values')
    ax[0].bar(ind + width/2, mean_3a, width, label =  'ADCC Values')
    ax[0].set_title("ADCC (Fig. 3A)")
    ax[0].set_xlabel("Glycans")
    ax[0].set_xticklabels(mixture, rotation=90)
    ax[0].set_xticks(ind)
    ax[0].legend()

    ax[1].bar(ind - width/2, infer_adcc_3b, width, label = 'Inferred ADCC Values')
    ax[1].bar(ind + width/2, mean_3b , width, label = 'ADCC Values')
    ax[1].set_title("ADCC (Fig. 3B)")
    ax[1].set_xticklabels(mixture, rotation=90)
    ax[1].set_xticks(ind)
    ax[1].legend()

    # Add subplot labels
    subplotLabel(ax)

    return f