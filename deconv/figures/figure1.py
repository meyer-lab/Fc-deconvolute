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

    ax, f = getSetup((10, 3), (1, 3))

    data_dekkers = load_dekkers()
    df = data_dekkers["profiling"]

    data = df.groupby(["index", "receptor"]).mean().reset_index()
    data2 = data.pivot(index="index", columns="receptor", values="binding")
    pca2 = PCA()
    data_new = pca2.fit_transform(data2)

    acc_variance = pca2.explained_variance_ratio_.copy()[:6]
    acc_variance = np.cumsum(acc_variance)

    ax[0].plot(range(1, acc_variance.size + 1), acc_variance)
    ax[0].set_ylabel("Explained Variance")
    ax[0].set_xlabel("Number of Components")
    ax[0].set_ylim([0, 1])

    ax[1].scatter(data_new[:, 0], data_new[:, 1])
    ax[1].set_title("Scores")
    ax[1].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[1].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[1]))

    mixtures = data_dekkers["mixtures"]
    for i in range(20):
        ax[1].annotate(mixtures[i], (data_new[i, 0], data_new[i, 1]))

    loadings = pd.DataFrame(pca2.components_.T[:, :2], columns=['PC1', 'PC2'], index=data2.columns)
    
    ax[2].set_title("Loadings")
    ax[2].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[2].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[1]))

    for i in range(12):
        ax[2].scatter(loadings.iloc[i, 0], loadings.iloc[i, 1], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')

    
    # Add subplot labels
    subplotLabel(ax)

    return f