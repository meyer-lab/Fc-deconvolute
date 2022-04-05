from deconv.imports import load_dekkers
from deconv.figures.common import subplotLabel, getSetup
import numpy as np
import seaborn as sns

def makeFigure():
    ax, f = getSetup((7, 5), (1,1))

    data_dekkers = load_dekkers()
    antiD = data_dekkers["antiD"]
    glycans = data_dekkers["glycans"]

    corr = np.corrcoef(np.transpose(antiD))
    mask = np.triu(np.ones_like(corr, dtype=bool))
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    ax[0] = sns.heatmap(corr, mask=mask, vmax=.3, center=0, square=True, cmap = cmap, linewidths=.5, cbar_kws={"shrink": .5}, xticklabels=glycans, yticklabels=glycans)
    
    # Add subplot labels
    subplotLabel(ax)

    return f
