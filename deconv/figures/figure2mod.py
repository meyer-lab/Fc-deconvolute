import numpy as np
import matplotlib.ticker as mticker
from .common import getSetup, subplotLabel
from ..emceeDeconv import getEmceeTrace
from ..imports import load_dekkers
import pandas as pd

def makeFigure():
    ax, f = getSetup((7.5, 7.5), (4, 3))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]

    data2 = data_dekkers["profiling"]

    trace = getEmceeTrace()

    activity = trace.posterior.activity[0]

    qqs = np.quantile(activity, (0.025, 0.5, 0.975), axis=0)
    median = qqs[1, :, :]
    p025 = median - qqs[0, :, :]
    p975 = qqs[2, :, :] - median

    g = ['G0N', 'G0', 'G0F', 'G0FN', 'G1', 'G1F', 'G1N', 'G1S', 'G1FS', 'G1FN', 'G1NS', 'G1FNS', 'G2', 'G2S', 'G2F', 'G2N', 'G2FN', 'G2NS', 'G2FS', 'G2FS2', 'G2S2', 'G2NS2', 'G2FNS', 'G2FNS2' ]

    median_df = pd.DataFrame(median)
    median_df.index = glycans
    median_df.columns = data2.columns
    median_df = median_df.reindex(g)

    p975_df = pd.DataFrame(p975)
    p975_df.index = glycans
    p975_df.columns = data2.columns
    p975_df = p975_df.reindex(g)

    p025_df = pd.DataFrame(p025)
    p025_df.index = glycans
    p025_df.columns = data2.columns
    p025_df = p025_df.reindex(g)

    for i in range(median.shape[1]):
        ax[i].errorbar(g, median_df.iloc[:, i], yerr=[p025_df.iloc[:, i], p975_df.iloc[:, i]], fmt='o', markersize = 5)
        ax[i].set_ylabel(median_df.columns[i], size=6)
        ax[i].xaxis.set_major_locator(mticker.FixedLocator(ax[i].get_xticks()))
        ax[i].set_xticklabels(g, rotation=90, size=6)
        _, _, _, y2 = ax[i].axis()
        ax[i].set_ylim([0, y2])

    subplotLabel(ax)

    return f