from sklearn.decomposition import PCA
from deconv.imports import load_dekkers
from deconv.figures.common import subplotLabel, getSetup
import pandas as pd
import numpy as np


def makeFigure():
    ax, f = getSetup((9, 3), (1, 3))

    data_dekkers = load_dekkers()

    mean_binding = data_dekkers["bindings"]
    mean_binding = [m.groupby(level=0).mean() for m in mean_binding]

    mean_3a = data_dekkers["meanADCC3a"]
    mean_3b = data_dekkers["meanADCC3b"]
    data = [mean_3a, mean_3b] + mean_binding

    # TODO: add compement activation

    data2 = np.transpose(np.array(data))
    pca2 = PCA()
    data_new = pca2.fit_transform(data2)

    ax[0].scatter(data_new[:, 0], data_new[:, 1])
    ax[0].set_title("Scores")
    ax[0].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[0].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[1]))

    l = ['adcc3a', 'adcc3b', 'bindingA', 'bindingB', 'bindingC', 'bindingD', 'bindingE', 'bindingF', 'bindingG', 'bindingH']
    loadings = pd.DataFrame(pca2.components_.T[:, :2], columns=['PC1', 'PC2'], index=l)
    ax[1].scatter(loadings.iloc[:, 0], loadings.iloc[:, 1])
    ax[1].set_title("Loadings")
    ax[1].set_xlabel("Component 1 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[0]))
    ax[1].set_ylabel("Component 2 ({ratio:.2f})".format(ratio=pca2.explained_variance_ratio_[1]))

    for i in range(10):
        ax[1].annotate(l[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    acc_variance = pca2.explained_variance_ratio_.copy()
    acc_variance = np.cumsum(acc_variance)

    ax[2].plot(range(1, acc_variance.size + 1), acc_variance)
    ax[2].set_ylabel("Explained Variance")
    ax[2].set_xlabel("Number of Components")

    # Add subplot labels
    subplotLabel(ax)

    return f
