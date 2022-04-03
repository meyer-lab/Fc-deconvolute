from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers


def makeFigure():
    #imports and formats
    trace = getEmceeTrace()

    activity_scores = (trace.posterior.activity_scores[0])
    activity_loadings = (trace.posterior.activity_loadings[0])
    median_scores = np.median(activity_scores, axis=0)
    median_loadings = np.transpose(np.median(activity_loadings, axis=0))

    # set up matrix w errors: SCORES
    lowErrS = np.subtract((np.percentile(activity_scores, 50, axis=0)), (np.percentile(activity_scores, 33, axis=0)))
    highErrS = np.subtract((np.percentile(activity_scores, 66, axis=0)), (np.percentile(activity_scores, 50, axis=0)))

    # set up matrix w errors: LOADINGS
    ac33 = (np.percentile(activity_loadings, 33, axis=0))
    ac66 = np.percentile(activity_loadings, 66, axis=0)
    ac33 = np.transpose(ac33)
    ac66 = np.transpose(ac66)
    lowErrL = np.subtract(median_loadings, ac33)
    highErrL = np.subtract(ac66, median_loadings)

    ax, f = getSetup((10, 7), (2, 3))

    ScoreMarkers = ['^','^','^', 'o','o', 'o', 'd','d', 'D', 's', 's', 'X', 'x', 'x', 'x', '+', '+','+', '1', '1', '2', '3', '3', '4']
    ScoreColor = ['orchid', 'darkblue', 'palegreen', 'orchid', 'darkblue', 'palegreen', 'darkblue', 'palegreen', 'palegreen', 'darkblue', 'palegreen', 'palegreen','orchid', 'darkblue', 'palegreen', 'orchid', 'darkblue', 'palegreen', 'darkblue', 'palegreen' , 'palegreen','darkblue', 'palegreen' , 'palegreen']

    ax[1].set_title("Activity Scores")
    ax[1].set_xlabel("Component 1")
    ax[1].set_ylabel("Component 2")
    for i in range(24):
        ax[1].scatter(median_scores[i, 0], median_scores[i, 1], marker = ScoreMarkers[i], color = ScoreColor[i], edgecolor='k')
    ax[1].errorbar(median_scores[:, 0], median_scores[:, 1], yerr=[lowErrS[:, 1], highErrS[:, 1]], xerr=[lowErrS[:, 0], highErrS[:, 0]], fmt = ',', color ='k', lw = .5)
    ax[1].set_ylim(bottom = 0)
    ax[1].set_xlim(left = 0)

    LoadingMarkers = ['o','o', '^', 'd', 's', 'v', 'v', 'v', 'o', 'o', 'o', 'o']
    LoadingColors = ['lightcoral', 'lightcoral', 'gold', 'gold', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise']

    ax[2].set_title("Activity Loadings")
    ax[2].set_xlabel("Component 1")
    ax[2].set_ylabel("Component 2")
    for i in range(12):
        ax[2].scatter(median_loadings[i, 0], median_loadings[i, 1], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k' )
    ax[2].errorbar(median_loadings[:, 0], median_loadings[:, 1], yerr=[lowErrL[:, 1], highErrL[:, 1]], xerr=[lowErrL[:, 0], highErrL[:, 0]], fmt = ',', color ='k', lw = .5)
    ax[2].set_ylim(bottom = 0)
    ax[2].set_xlim(left = 0)


    ax[4].set_title("Activity Scores")
    ax[4].set_xlabel("Component 1")
    ax[4].set_ylabel("Component 3")
    for i in range(24):
        ax[4].scatter(median_scores[i, 0], median_scores[i, 2], marker = ScoreMarkers[i], color = ScoreColor[i],edgecolor='k')
    ax[4].errorbar(median_scores[:, 0], median_scores[:, 2], yerr=[lowErrS[:, 2], highErrS[:, 2]], xerr=[lowErrS[:, 0], highErrS[:, 0]],fmt = ',', color ='k', lw = .5)
    ax[4].set_ylim(bottom = 0)
    ax[4].set_xlim(left = 0)

    ax[5].set_title("Activity Loadings")
    ax[5].set_xlabel("Component 1")
    ax[5].set_ylabel("Component 3")
    for i in range(12):
        ax[5].scatter(median_loadings[i, 0], median_loadings[i, 2], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')
    ax[5].errorbar(median_loadings[:, 0], median_loadings[:, 2], yerr=[lowErrL[:, 2], highErrL[:, 2]], xerr=[lowErrL[:, 0], highErrL[:, 0]],fmt = ',', color ='k', lw = .5)
    ax[5].set_ylim(bottom = 0)
    ax[5].set_xlim(left = 0)

    # Add subplot labels
    subplotLabel(ax)

    return f
