import numpy as np
from .common import getSetup, subplotLabel
from ..emceeDeconv import getEmceeTrace
from ..imports import load_dekkers


def makeFigure():
    ax, f = getSetup((7.5, 7.5), (2, 2))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]

    activity = getEmceeTrace()

    qqs = np.quantile(activity, (0.025, 0.5, 0.975), axis=0)
    median = qqs[1, :, :]
    p33 = median - qqs[0, :, :]
    p66 = qqs[2, :, :] - median

    colors = ["dodgerblue"] * 12 + ["salmon"] * 12

    i = 0
    ax[i].errorbar(median[:, 0], median[:, 1], xerr=[p33[:, 0], p66[:, 0]], yerr=[p33[:, 1], p66[:, 1]], fmt="none", ecolor=colors)
    ax[i].scatter(median[:, 0], median[:, 1], marker="o", color=colors)

    # zip joins x and y coordinates in pairs
    for j in range(median.shape[0]):
        ax[i].annotate(glycans[j], # this is the text
                       (median[j, 0], median[j, 1]), # these are the coordinates to position the label
                       textcoords="offset points", # how to position the text
                       xytext=(0, 10), # distance from text to points (x,y)
                       ha='center') # horizontal alignment can be left, right or center

    subplotLabel(ax)

    return f
