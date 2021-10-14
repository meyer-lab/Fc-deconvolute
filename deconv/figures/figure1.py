from sklearn.decomposition import NMF
from deconv.imports import load_figures, load_bindingData
from deconv.figures.common import subplotLabel, getSetup
import pandas as pd
import numpy as np


def makeFigure():
    ax, f = getSetup((9, 3), (1, 3))

    adcc3a, adcc3b = load_figures()
    mean_binding = list(load_bindingData())
    mean_binding = [m.groupby(level=0).mean() for m in mean_binding]

    mean_3a = (adcc3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc3b.groupby(level=0).sum()) / 4
    data = [mean_3a, mean_3b] + mean_binding

    data2 = np.transpose(np.array(data))
    pca2 = NMF(n_components=2, max_iter=2000, tol=1e-9, init="nndsvda")
    data_new = pca2.fit_transform(data2)

    ax[0].scatter(data_new[:, 0], data_new[:, 1])
    ax[0].set_title("Scores")
    ax[0].set_xlabel("Component 1")
    ax[0].set_ylabel("Component 2")

    l = ['adcc3a', 'adcc3b', 'bindingA', 'bindingB', 'bindingC', 'bindingD', 'bindingE', 'bindingF', 'bindingG', 'bindingH']
    loadings = pd.DataFrame(pca2.components_.T[:, :2], columns=['PC1', 'PC2'], index=l)
    ax[1].scatter(loadings.iloc[:, 0], loadings.iloc[:, 1])
    ax[1].set_title("Loadings")
    ax[1].set_xlabel("Component 1")
    ax[1].set_ylabel("Component 2")

    for i in range(10):
        ax[1].annotate(l[i], (loadings.iloc[i, 0], loadings.iloc[i, 1]))

    error = np.zeros(5)
    for ii in range(5):
        pca2 = NMF(n_components=ii+1, max_iter=2000, tol=1e-9, init="nndsvda")
        pca2.fit(data2)

        error[ii] = pca2.reconstruction_err_

    ax[2].plot(range(1, 6), error)
    ax[2].set_ylabel("Error")
    ax[2].set_xlabel("Number of Components")
    ax[2].set_ylim(bottom=0.0)

    # Add subplot labels
    subplotLabel(ax)

    return f
