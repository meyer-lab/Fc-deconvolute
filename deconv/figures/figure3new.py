from sklearn.decomposition import PCA
from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers

def makeFigure():
    ax, f = getSetup((12, 6), (1, 2))
    trace = getEmceeTrace()
    data_dekkers = load_dekkers()

    activity_scores = np.median(trace.posterior.activity_scores[0], axis =0)
    activity_loadings = np.median(trace.posterior.activity_loadings[0], axis = 0)

    pca = PCA()
    pca2 = PCA()
    scores_pca = pca.fit_transform(activity_scores)
    loadings_pca = pca2.fit_transform(activity_loadings)

    ax, f = getSetup((9, 4), (1, 2))

    ax[0].scatter(scores_pca[:, 0], scores_pca[:, 1])
    ax[0].set_title("Activity Scores")
    ax[0].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pca.explained_variance_ratio_[0]))
    ax[0].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pca.explained_variance_ratio_[1]))
    glycans = data_dekkers["glycans"]

    scores = pd.DataFrame(scores_pca)
    for i in range(24):
        ax[0].annotate(glycans[i], (scores.iloc[i, 0], scores.iloc[i, 1]))

    ax[1].scatter(loadings_pca[:, 0], loadings_pca[:, 1])
    ax[1].set_title("Activity Loadings")
    ax[1].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[1].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[1]))
    labels = ["FcRI", 'FcRII', 'FcRIII']
    loadings = pd.DataFrame(loadings_pca)

    for i in range(3):
        ax[1].annotate(labels[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    # Add subplot labels
    subplotLabel(ax)
    return f