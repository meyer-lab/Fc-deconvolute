from deconv.imports import load_dekkers
from deconv.figures.common import subplotLabel, getSetup
from deconv.pca import pca
import pandas as pd
import numpy as np


def makeFigure():
    # set loadings marker color and shape
    LoadingMarkers = ["o"] * 6 + ["v", "v", "v", "s", "^", "d"]
    LoadingColors = ["lightcoral"] * 2 + ["mediumturquoise"] * 8 + ["gold"] * 2

    ax, f = getSetup((14, 7), (2, 4))

    data_dekkers = load_dekkers()
    data2 = data_dekkers["profiling"]

    mixture = data_dekkers["antiD"]
    data_new, components_, expVar = pca(mixture)
    print(expVar)

    acc_variance = np.cumsum(expVar)
    xx = np.arange(1, acc_variance.size + 1)
    acc_variance = np.insert(acc_variance, 0, 0)
    xx = np.insert(xx, 0, 0)

    ax[0].plot(xx, acc_variance)
    ax[0].plot(xx, xx / mixture.shape[0])
    ax[0].set_ylabel("Explained Variance")
    ax[0].set_xlabel("Number of Components")
    ax[0].set_title("Mixture PCA")
    ax[0].set_ylim([0, 1])

    data_new, components_, expVar = pca(data2)

    acc_variance = np.cumsum(expVar)[:6]

    ax[1].plot(range(1, acc_variance.size + 1), acc_variance)
    ax[1].set_ylabel("Explained Variance")
    ax[1].set_xlabel("Number of Components")
    ax[1].set_ylim([0, 1])

    ax[2].scatter(data_new[:, 0], data_new[:, 1])
    ax[2].set_title("Scores")
    ax[2].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=expVar[0]))
    ax[2].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=expVar[1]))

    mixtures = data_dekkers["mixtures"]
    for i in range(20):
        ax[2].annotate(mixtures[i], (data_new[i, 0], data_new[i, 1]))

    loadings = pd.DataFrame(components_.T[:, :3], columns=["PC1", "PC2", "PC3"], index=data2.columns)

    ax[3].set_title("Loadings")
    ax[3].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=expVar[0]))
    ax[3].set_ylabel("Component 2 ({ratio:.0%})".format(ratio=expVar[1]))

    for i in range(12):
        ax[3].scatter(loadings.iloc[i, 0], loadings.iloc[i, 1], marker=LoadingMarkers[i], color=LoadingColors[i], edgecolor="k")

    ax[6].scatter(data_new[:, 0], data_new[:, 2])
    ax[6].set_title("Scores")
    ax[6].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=expVar[0]))
    ax[6].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=expVar[2]))

    mixtures = data_dekkers["mixtures"]
    for i in range(20):
        ax[6].annotate(mixtures[i], (data_new[i, 0], data_new[i, 2]))

    ax[7].set_title("Loadings")
    ax[7].set_xlabel("Component 1 ({ratio:.0%})".format(ratio=expVar[0]))
    ax[7].set_ylabel("Component 3 ({ratio:.0%})".format(ratio=expVar[2]))

    for i in range(12):
        ax[7].scatter(loadings.iloc[i, 0], loadings.iloc[i, 2], marker=LoadingMarkers[i], color=LoadingColors[i], edgecolor="k")

    # Add subplot labels
    subplotLabel(ax)

    return f
