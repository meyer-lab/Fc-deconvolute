from sklearn.decomposition import PCA
from deconv.imports import load_dekkers
from deconv.figures.common import subplotLabel, getSetup
import pandas as pd
import numpy as np


def makeFigure():
    # set loadings marker color and shape
    LoadingMarkers = ['o','o', 'o', 'o', 'o', 'o', 'v', 
        'v', 'v', 's', '^', 'd']
    LoadingColors =['lightcoral', 'lightcoral',  
        'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 
        'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 
        'mediumturquoise', 'mediumturquoise','gold', 'gold']

    ax, f = getSetup((14, 7), (2, 4))

    data_dekkers = load_dekkers()
    data2 = data_dekkers["profiling"]

    pca2 = PCA()
    data_new = pca2.fit_transform(data2)

    acc_variance = pca2.explained_variance_ratio_.copy()[:6]
    acc_variance = np.cumsum(acc_variance)

    ax[1].plot(range(1, acc_variance.size + 1), acc_variance)
    ax[1].set_ylabel("Explained Variance")
    ax[1].set_xlabel("Number of Components")
    ax[1].set_ylim([0, 1])

    ax[2].scatter(data_new[:, 0], data_new[:, 1])
    ax[2].set_title("Scores")
    ax[2].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[2].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[1]))

    mixtures = data_dekkers["mixtures"]
    for i in range(20):
        ax[2].annotate(mixtures[i], (data_new[i, 0], data_new[i, 1]))

    loadings = pd.DataFrame(pca2.components_.T[:, :2], columns=['PC1', 'PC2'], index=data2.columns)
    
    ax[3].set_title("Loadings")
    ax[3].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[3].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[1]))

    for i in range(12):
        ax[3].scatter(loadings.iloc[i, 0], loadings.iloc[i, 1], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')

    ax[6].scatter(data_new[:, 0], data_new[:, 2])
    ax[6].set_title("Scores")
    ax[6].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[6].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[2]))

    mixtures = data_dekkers["mixtures"]
    for i in range(20):
        ax[6].annotate(mixtures[i], (data_new[i, 0], data_new[i, 2]))

    loadings = pd.DataFrame(pca2.components_.T[:, [0, 2]], columns=['PC1', 'PC3'], index=data2.columns)
    ax[7].set_title("Loadings")
    ax[7].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[7].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[2]))

    for i in range(12):
        ax[7].scatter(loadings.iloc[i, 0], loadings.iloc[i, 1], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')

    # Add subplot labels
    subplotLabel(ax)

    return f