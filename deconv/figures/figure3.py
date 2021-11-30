from sklearn.decomposition import PCA
from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers

def makeFigure():
    #imports and formats
    ax, f = getSetup((12, 6), (1, 2))
    trace = getEmceeTrace()
    data_dekkers = load_dekkers()
    pca = PCA()
    pca2 = PCA()

    activity_scores = (trace.posterior.activity_scores[0])
    activity_loadings = (trace.posterior.activity_loadings[0])
    median_scores = np.median(activity_scores, axis = 0)
    median_loadings = np.median(activity_loadings, axis = 0)

    scorespca = pca.fit(median_scores)
    loadingspca =pca2.fit(median_loadings)
    pc1_var = pca.explained_variance_ratio_[0]
    pc2_var = pca.explained_variance_ratio_[1]
    p1var = pca2.explained_variance_ratio_[0]
    p2var = pca2.explained_variance_ratio_[1]

    #set up matrix w errors: SCORES
    p33_s = np.zeros((24,3))
    p66_s = np.zeros((24,3))
    for i in range(3):
        p33_s[:, i] = np.percentile(median_scores[:,i], 33)
        p66_s[:,i] =np.percentile(median_scores[:,i], 66)
    

    #set up matrix w errors: LOADINGS
    median_loadings = np.transpose(median_loadings)
    p33_l = np.zeros((12,3))
    p66_l = np.zeros((12,3))
    for i in range(3):
        p33_l[:, i] = np.percentile(median_loadings[:,i], 33)
        p66_l[:,i] =np.percentile(median_loadings[:,i], 66)

    ax, f = getSetup((9, 4), (1, 2))

    ax[0].set_title("Activity Scores")
    ax[0].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pc1_var))
    ax[0].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pc2_var))
    ax[0].errorbar(median_scores[:, 0], median_scores[:, 1], yerr = [p33_s[:,1], p66_s[:,1]], xerr= [p33_s[:,0], p66_s[:,0]], fmt='o')
    glycans = data_dekkers["glycans"]

    scores = pd.DataFrame(median_scores)
    for i in range(24):
        ax[0].annotate(glycans[i], (scores.iloc[i, 0], scores.iloc[i, 1]))

    ax[1].set_title("Activity Loadings")
    ax[1].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=p1var))
    ax[1].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=p2var))
    ax[1].errorbar(median_loadings[:, 0], median_loadings[:, 1], yerr = [p33_l[:,1], p66_l[:,1]], xerr = [p33_l[:,0], p66_l[:,0]], fmt = 'o')
    labels = ['ADCC FcγRIIIA158F/F', 'ADCC FcγRIIIA158V/V', 'Complement Activation C1q', 'Complement Activation C4', 'Binding FcγRIa', 'Binding FcγRIIa 131H', 'Binding FcγRIIa 131R', 'Binding FcγRIIb/c', 'Binding FcγRIIIa 158F', 'Binding FcγRIIIa 158V', 'Binding Fc-FcγRIIIb NA1', 'Binding Fc-FcγRIIIb NA2']
    loadings = pd.DataFrame(median_loadings)

    for i in range(12):
        ax[1].annotate(labels[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    # Add subplot labels
    subplotLabel(ax)
    
    return f