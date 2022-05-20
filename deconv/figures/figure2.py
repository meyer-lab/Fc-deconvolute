import numpy as np
import matplotlib.ticker as mticker
from .common import getSetup, subplotLabel
from ..emceeDeconv import getEmceeTrace
from ..imports import load_dekkers


def makeFigure():
    ax, f = getSetup((7.5, 7.5), (4, 3))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]

    data = data_dekkers["profiling"]
    data2 = data.pivot(index="index", columns="receptor", values="binding")

    trace = getEmceeTrace()

    activity = trace.posterior.activity[0]

    qqs = np.quantile(activity, (0.025, 0.5, 0.975), axis=0)
    median = qqs[1, :, :]
    p33 = median - qqs[0, :, :]
    p66 = qqs[2, :, :] - median

    for i in range(median.shape[1]):
        ax[i].errorbar(glycans, median[:, i], yerr=[p33[:, i], p66[:, i]], fmt='o', markersize = 5)
        ax[i].set_ylabel(data2.columns[i], size=6)
        ax[i].xaxis.set_major_locator(mticker.FixedLocator(ax[i].get_xticks()))
        ax[i].set_xticklabels(glycans, rotation=90, size=6)
        _, _, _, y2 = ax[i].axis()
        ax[i].set_ylim([0, y2])

    subplotLabel(ax)

    return f
