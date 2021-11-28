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

    activity_scores = (trace.posterior.activity_scores[0])
    activity_loadings = (trace.posterior.activity_loadings[0])
    median_scores = np.median(activity_scores, axis = 0)
    median_loadings = np.median(activity_loadings, axis = 0)
    
    pca = PCA()
    pca2 = PCA()
    scores_pca = pca.fit_transform(median_scores)
    loadings_pca = pca2.fit_transform(median_loadings)

    #set up matrix w errors: SCORES
    p33_s = np.zeros((24,3))
    p66_s = np.zeros((24,3))
    for i in range(3):
        p33_s[:, i] = np.percentile(median_scores[:,i], 33)
        p66_s[:,i] =np.percentile(median_scores[:,i], 66)

    err_pca_33s = pca.transform(p33_s)
    err_pca_66s = pca.transform(p66_s)

    #set up matrix w errors: LOADINGS
    p33_l = np.zeros((3,12))
    p66_l = np.zeros((3,12))
    for i in range(12):
        p33_l[:, i] = np.percentile(median_loadings[:,i], 33)
        p66_l[:,i] =np.percentile(median_loadings[:,i], 66)

    err_pca_33l = pca2.transform(p33_l)
    err_pca_66l = pca2.transform(p66_l)


    ax, f = getSetup((9, 4), (1, 2))

    #ax[0].scatter(scores_pca[:, 0], scores_pca[:, 1])
    ax[0].set_title("Activity Scores")
    ax[0].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pca.explained_variance_ratio_[0]))
    ax[0].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pca.explained_variance_ratio_[1]))
    ax[0].errorbar(scores_pca[:, 0], scores_pca[:, 1], yerr=[err_pca_33s[:,1], err_pca_66s[:,1]], xerr = [err_pca_33s[:,0], err_pca_66s[:,0]], fmt='o')
    glycans = data_dekkers["glycans"]

    scores = pd.DataFrame(scores_pca)
    for i in range(24):
        ax[0].annotate(glycans[i], (scores.iloc[i, 0], scores.iloc[i, 1]))

    #ax[1].scatter(loadings_pca[:, 0], loadings_pca[:, 1])
    ax[1].set_title("Activity Loadings")
    ax[1].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[1].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[1]))
    ax[1].errorbar(loadings_pca[:, 0], loadings_pca[:, 1], yerr=[err_pca_33l[:,1], err_pca_66l[:,1]], xerr = [err_pca_33l[:,0], err_pca_66l[:,0]], fmt='o')
    labels = ["FcRI", 'FcRII', 'FcRIII']
    loadings = pd.DataFrame(loadings_pca)

    for i in range(3):
        ax[1].annotate(labels[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    # Add subplot labels
    subplotLabel(ax)
    
    return f