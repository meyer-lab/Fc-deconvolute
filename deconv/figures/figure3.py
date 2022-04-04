from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers


def makeFigure():
    #imports and formats
    trace = getEmceeTrace()
    data_dekkers = load_dekkers()

    df = data_dekkers["profiling"]
    data = df.groupby(["index", "receptor"]).mean().reset_index()
    data2 = data.pivot(index="index", columns="receptor", values="binding")

    activity_scores = (trace.posterior.activity_scores[0])
    activity_loadings = (trace.posterior.activity_loadings[0])
    median_scores = np.median(activity_scores, axis=0)
    median_loadings = np.transpose(np.median(activity_loadings, axis=0))

    # set up matrix w errors: SCORES
    lowErrS = np.subtract((np.percentile(activity_scores, 50, axis=0)), (np.percentile(activity_scores, 33, axis=0)))
    highErrS = np.subtract((np.percentile(activity_scores, 66, axis=0)), (np.percentile(activity_scores, 50, axis=0)))

    # set up matrix w errors: LOADINGS
    ac33 = np.percentile(activity_loadings, 33, axis=0)
    ac66 = np.percentile(activity_loadings, 66, axis=0)
    ac33 = np.transpose(ac33)
    ac66 = np.transpose(ac66)
    lowErrL = np.subtract(median_loadings, ac33)
    highErrL = np.subtract(ac66, median_loadings)

    ax, f = getSetup((7, 7), (2, 2))

    ax[0].set_title("Activity Scores")
    ax[0].set_xlabel("Component 1")
    ax[0].set_ylabel("Component 2")
    ax[0].set_ylim([0, 1.5])
    ax[0].set_xlim([0, 7])
    ax[0].errorbar(median_scores[:, 0], median_scores[:, 1], yerr=[lowErrS[:, 1], highErrS[:, 1]], xerr=[lowErrS[:, 0], highErrS[:, 0]], fmt='o')
    glycans = data_dekkers["glycans"]

    scores = pd.DataFrame(median_scores)
    for i in range(24):
        ax[0].annotate(glycans[i], (scores.iloc[i, 0], scores.iloc[i, 1]))

    ax[1].set_title("Activity Loadings")
    ax[1].set_xlabel("Component 1")
    ax[1].set_ylabel("Component 2")
    ax[1].set_ylim([0, 1])
    ax[1].set_xlim([0, 1])
    ax[1].errorbar(median_loadings[:, 0], median_loadings[:, 1], yerr=[lowErrL[:, 1], highErrL[:, 1]], xerr=[lowErrL[:, 0], highErrL[:, 0]], fmt='o')
    labels = data2.columns
    loadings = pd.DataFrame(median_loadings)

    for i in range(12):
        ax[1].annotate(labels[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    ax[2].set_title("Activity Scores")
    ax[2].set_xlabel("Component 1")
    ax[2].set_ylabel("Component 3")
    ax[2].set_ylim([0, 5])
    ax[2].set_xlim([0, 7])
    ax[2].errorbar(median_scores[:, 0], median_scores[:, 2], yerr=[lowErrS[:, 2], highErrS[:, 2]], xerr=[lowErrS[:, 0], highErrS[:, 0]], fmt='o')

    scores = pd.DataFrame(median_scores)
    for i in range(24):
        ax[2].annotate(glycans[i], (scores.iloc[i, 0], scores.iloc[i, 2]))

    ax[3].set_title("Activity Loadings")
    ax[3].set_xlabel("Component 1")
    ax[3].set_ylabel("Component 3")
    ax[3].set_ylim([0, 0.3])
    ax[3].set_xlim([0, 1])
    ax[3].errorbar(median_loadings[:, 0], median_loadings[:, 2], yerr=[lowErrL[:, 2], highErrL[:, 2]], xerr=[lowErrL[:, 0], highErrL[:, 0]], fmt='o')

    for i in range(12):
        ax[3].annotate(labels[i], (loadings.iloc[i, 0], loadings.iloc[i, 2]))

    # Add subplot labels
    subplotLabel(ax)

    return f
