from deconv.imports import load_dekkers
from deconv.figures.common import subplotLabel, getSetup
import numpy as np
import seaborn as sns

def makeFigure():
    ax, f = getSetup((7, 5), (1,1))

    data_dekkers = load_dekkers()
    antiD = data_dekkers["antiD"]
    glycans = data_dekkers["glycans"]
    mixtures = data_dekkers["mixtures"]

    ax[0] = sns.heatmap(antiD, linewidth=0.5, xticklabels=glycans, yticklabels=mixtures)
    ax[0].collections[0].colorbar.set_label("Relative Abundance of Glycopeptides")

     # Add subplot labels
    subplotLabel(ax)

    return f