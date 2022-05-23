from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


def makeFigure():
    #imports and formats
    activity = getEmceeTrace()

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

    plotPCA(ax[0], median_scores, lowErrS, highErrS, ScoreMarkers, ScoreColor, 1)
    ax[0].set_title("Activity Scores")

    legend1 = ax[0].legend(handles= legend_elements1, bbox_to_anchor=(0., 1.25, 1., .25), loc='lower left',
                      ncol=1, borderaxespad=0., prop= {'size': 11}, handlelength = 2)
    ax[0].add_artist(legend1)

    plotPCA(ax[1], median_loadings, lowErrL, highErrL, LoadingMarkers, LoadingColors, 1)
    ax[1].set_title("Activity Loadings")

    legend2 = ax[1].legend(handles= legend_elements2, bbox_to_anchor=(0., 1.25, .7, .25), loc='lower right',
                      ncol=1, borderaxespad=0., prop= {'size': 12}, handleheight = 2)
    ax[1].add_artist(legend2)

    ax[2].set_title("Activity Scores")
    plotPCA(ax[2], median_scores, lowErrS, highErrS, ScoreMarkers, ScoreColor, 2)

    ax[3].set_title("Activity Loadings")
    plotPCA(ax[3], median_loadings, lowErrL, highErrL, LoadingMarkers, LoadingColors, 2)

    # Add subplot labels
    subplotLabel(ax)

    return f


def plotPCA(ax, median, lowError, highError, ScoreMarkers, ScoreColor, compTwo: int):
    """Plot the scores/loadings."""
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component " + str(compTwo + 1))
    ax.errorbar(median[:, 0], median[:, compTwo],
         yerr=[lowError[:, compTwo], highError[:, compTwo]], xerr=[lowError[:, 0], 
         highError[:, 0]],  fmt = ',', color ='k', lw = .25)

    for i in range(median.shape[0]):
        ax.scatter(median[i, 0], median[i, compTwo], marker = ScoreMarkers[i], color = ScoreColor[i])