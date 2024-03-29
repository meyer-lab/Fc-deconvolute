import numpy as np
import matplotlib.ticker as mticker
from .common import getSetup, subplotLabel
from ..emceeDeconv import getEmceeTrace
from ..imports import load_dekkers


def makeFigure():
    ax, f = getSetup((7.5, 7.5), (4, 3))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]

    data2 = data_dekkers["profiling"]

    activity = getEmceeTrace()

    qqs = np.quantile(activity, (0.025, 0.5, 0.975), axis=0)
    median = qqs[1, :, :]
    p33 = median - qqs[0, :, :]
    p66 = qqs[2, :, :] - median

    # Difference with fucosylation
    fucDiff = activity[:, 12:, :] / activity[:, :12, :]
    fucQ = np.quantile(fucDiff, (0.025, 0.5, 0.975), axis=0)

    colors = ["dodgerblue"] * 12 + ["salmon"] * 12

    for i in range(median.shape[1]):
        ax[i].errorbar(glycans, median[:, i], yerr=[p33[:, i], p66[:, i]], fmt="none", ecolor=colors)
        ax[i].scatter(glycans, median[:, i], marker="o", color=colors)
        ax[i].set_ylabel(data2.columns[i], size=6)
        ax[i].xaxis.set_major_locator(mticker.FixedLocator(ax[i].get_xticks()))
        ax[i].set_xticklabels(glycans, rotation=90, size=6)
        if i < 6:
            ax[i].set_ylim([0, 8])
        else:
            ax[i].set_ylim([0, 4])

    subplotLabel(ax)

    return f
