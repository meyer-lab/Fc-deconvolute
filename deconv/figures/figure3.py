from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers
from sklearn.decomposition import PCA
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


def makeFigure():
    #imports and formats
    trace = getEmceeTrace()
    data_dekkers = load_dekkers()

    df = data_dekkers["profiling"]
    data = df.groupby(["index", "receptor"]).mean().reset_index()
    data2 = data.pivot(index="index", columns="receptor", values="binding")

    activity = trace.posterior.activity[0]

    activity_scores = []
    activity_loadings = []

    pca2 = PCA()

    for ii in range(activity.shape[0]):
        activity_scores.append(pca2.fit_transform(np.squeeze(activity[ii, :, :])))
        activity_loadings.append(pca2.components_)

    median_scores = np.median(activity_scores, axis=0)
    median_loadings = np.median(activity_loadings, axis=0).T

    # set up matrix w errors: SCORES
    lowErrS = np.subtract(median_scores, np.percentile(activity_scores, 33, axis=0))
    highErrS = np.subtract(np.percentile(activity_scores, 66, axis=0), median_scores)

    # set up matrix w errors: LOADINGS
    ac33 = np.percentile(activity_loadings, 33, axis=0)
    ac66 = np.percentile(activity_loadings, 66, axis=0)
    lowErrL = np.subtract(median_loadings, ac33.T)
    highErrL = np.subtract(ac66.T, median_loadings)

    # set up color and shape markers
    ScoreMarkers = ['^','^','^', 'o','o', 'o', 'd','d', 
        'D', 's', 's', 'X', 'x', 'x', 'x', '+', '+','+', '1', 
        '1', '2', '|', '|', '_']
    ScoreColor = ['crimson', 'orchid', 'blue', 'crimson', 
         'orchid', 'blue', 'orchid', 'blue', 'blue',
         'orchid', 'blue', 'blue','crimson', 'orchid', 
         'blue', 'crimson', 'orchid', 'blue', 'orchid', 
         'blue' , 'blue','orchid', 'blue' , 'blue']
    LoadingMarkers = ['o','o', 'o', 'o', 'o', 'o', 'v', 
        'v', 'v', 's', '^', 'd']
    LoadingColors =['lightcoral', 'lightcoral',  
        'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 
        'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 
        'mediumturquoise', 'mediumturquoise','gold', 'gold']

    #set up legend elements
    legend_elements1 = [mlines.Line2D([], [], marker='^', markeredgecolor='k', label='Complement C1q',markersize=12, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='d', markeredgecolor='k', label='Complement C4',markersize=12, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='s', markeredgecolor='k', label='FcγRIa',markersize=12, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='v', markeredgecolor='k', label='FcγRII',markersize=12, ls = 'none', markerfacecolor = 'none'),
                   mlines.Line2D([], [], marker='o', markeredgecolor='k', label='FcγRIII',markersize=12, ls = 'none', markerfacecolor=  'none')]


    legend_elements2 = [mpatches.Patch(facecolor = 'lightcoral', label = 'ADCC'),
                   mpatches.Patch(facecolor = 'gold', label = 'Complement Activation'),
                   mpatches.Patch(facecolor = 'mediumturquoise', label = 'Binding')]


    ax, f = getSetup((7, 7), (2, 2))

    ax[0].set_title("Activity Scores")
    ax[0].set_xlabel("Component 1")
    ax[0].set_ylabel("Component 2")
    ax[0].errorbar(median_scores[:, 0], median_scores[:, 1],
         yerr=[lowErrS[:, 1], highErrS[:, 1]], xerr=[lowErrS[:, 0], 
         highErrS[:, 0]],  fmt = ',', color ='k', lw = .25)
    for i in range(24):
        ax[0].scatter(median_scores[i, 0], median_scores[i, 1], marker = ScoreMarkers[i], color = ScoreColor[i], edgecolor='k')

    legend1 = ax[0].legend(handles= legend_elements1, bbox_to_anchor=(0., 1.25, 1., .25), loc='lower left',
                      ncol=1, borderaxespad=0., prop= {'size': 11}, handlelength = 2)

    ax[0].add_artist(legend1)

    ax[1].set_title("Activity Loadings")
    ax[1].set_xlabel("Component 1")
    ax[1].set_ylabel("Component 2")
    ax[1].errorbar(median_loadings[:, 0], median_loadings[:, 1],
        yerr=[lowErrL[:, 1], highErrL[:, 1]], xerr=[lowErrL[:, 0], 
        highErrL[:, 0]], fmt=',', color = 'k', lw = .25)

    legend2 = ax[1].legend(handles= legend_elements2, bbox_to_anchor=(0., 1.25, .7, .25), loc='lower right',
                      ncol=1, borderaxespad=0., prop= {'size': 12}, handleheight = 2)
    ax[1].add_artist(legend2)

    for i in range(12):
        ax[1].scatter(median_loadings[i, 0], median_loadings[i, 1], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')

    ax[2].set_title("Activity Scores")
    ax[2].set_xlabel("Component 1")
    ax[2].set_ylabel("Component 3")
    ax[2].errorbar(median_scores[:, 0], median_scores[:, 2], 
        yerr=[lowErrS[:, 2], highErrS[:, 2]], xerr=[lowErrS[:, 0], 
        highErrS[:, 0]], fmt=',', color = 'k', lw = .25)

    for i in range(24):
        ax[2].scatter(median_scores[i, 0], median_scores[i, 2], marker = ScoreMarkers[i], color = ScoreColor[i], edgecolor='k')

    ax[3].set_title("Activity Loadings")
    ax[3].set_xlabel("Component 1")
    ax[3].set_ylabel("Component 3")
    ax[3].errorbar(median_loadings[:, 0], median_loadings[:, 2], 
        yerr=[lowErrL[:, 2], highErrL[:, 2]], xerr=[lowErrL[:, 0], 
        highErrL[:, 0]], fmt=',', color = 'k', lw = .25)

    for i in range(12):
        ax[3].scatter(median_loadings[i, 0], median_loadings[i, 2], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')

    # Add subplot labels
    subplotLabel(ax)

    return f