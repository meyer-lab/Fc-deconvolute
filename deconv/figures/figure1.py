from sklearn.decomposition import PCA
from deconv.imports import load_dekkers
from deconv.figures.common import subplotLabel, getSetup
import pandas as pd
import numpy as np


def makeFigure():
    ax, f = getSetup((15, 8), (2, 4))

    data_dekkers = load_dekkers()
    df = data_dekkers["profiling"]

    data = df.groupby(["index", "receptor"]).mean().reset_index()
    data2 = data.pivot(index="index", columns="receptor", values="binding")
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

    l = ['ADCC FcγRIIIA158F/F', 'ADCC FcγRIIIA158V/V', 'Comp Act C1q', 'Comp Act C4', 'Binding FcγRIa', 'Binding FcγRIIa 131H', 'Binding FcγRIIa 131R', 'Binding FcγRIIb/c', 'Binding FcγRIIIa 158F', 'Binding FcγRIIIa 158V', 'Binding Fc-FcγRIIIb NA1', 'Binding Fc-FcγRIIIb NA2']
    loadings = pd.DataFrame(pca2.components_.T[:, :2], columns=['PC1', 'PC2'], index=l)
    LoadingMarkers = ['o','o', '^', 'd', 's', 'v', 'v', 'v', 'o', 'o', 'o', 'o']
    LoadingColors =['lightcoral', 'lightcoral', 'gold', 'gold', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise', 'mediumturquoise']
    for i in range(12):
        ax[3].scatter(loadings.iloc[i, 0], loadings.iloc[i, 1], marker = LoadingMarkers[i], color = LoadingColors[i], edgecolor='k')    
        
    loadings = pd.DataFrame(pca2.components_.T[:, :2], columns=['PC1', 'PC2'], index=data2.columns)
    ax[3].scatter(loadings.iloc[:, 0], loadings.iloc[:, 1])
    ax[3].set_title("Loadings")
    ax[3].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[3].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[1]))

    for i in range(12):
        ax[3].annotate(data2.columns[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    ax[5].scatter(data_new[:, 0], data_new[:, 2])
    ax[5].set_title("Scores")
    ax[5].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[5].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[2]))

    mixtures = data_dekkers["mixtures"]
    for i in range(20):
        ax[5].annotate(mixtures[i], (data_new[i,0], data_new[i,2]))

    l = ['ADCC FcγRIIIA158F/F', 'ADCC FcγRIIIA158V/V', 'Comp Act C1q', 'Comp Act C4', 'Binding FcγRIa', 'Binding FcγRIIa 131H', 'Binding FcγRIIa 131R', 'Binding FcγRIIb/c', 'Binding FcγRIIIa 158F', 'Binding FcγRIIIa 158V', 'Binding Fc-FcγRIIIb NA1', 'Binding Fc-FcγRIIIb NA2']
    loadings = pd.DataFrame(pca2.components_.T[:, [0,2]], columns=['PC1', 'PC3'], index=l)
    for i in range(12):
        ax[6].scatter(loadings.iloc[i, 0], loadings.iloc[i, 1],color = LoadingColors[i], marker = LoadingMarkers[i], edgecolor='k')    
    ax[6].set_title("Loadings")
    ax[6].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[6].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[2]))

    ax[4].annotate(mixtures[i], (data_new[i, 0], data_new[i, 2]))

    loadings = pd.DataFrame(pca2.components_.T[:, [0, 2]], columns=['PC1', 'PC3'], index=data2.columns)
    ax[5].scatter(loadings.iloc[:, 0], loadings.iloc[:, 1])
    ax[5].set_title("Loadings")
    ax[5].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[5].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=pca2.explained_variance_ratio_[2]))

    for i in range(12):
        ax[5].annotate(data2.columns[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    # Add subplot labels
    subplotLabel(ax)

    return f
