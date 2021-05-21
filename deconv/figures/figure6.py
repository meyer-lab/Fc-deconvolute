from sklearn.decomposition import PCA
from deconv.imports import load_tables, load_figures, load_bindingData
import matplotlib.pyplot as plt
from deconv.figures.common import subplotLabel, getSetup
import pandas as pd
import numpy as np

def makeFigure():
    ax, f = getSetup((9, 3), (1, 3))

    pca2 = PCA(n_components=5)

    adcc3a, adcc3b = load_figures()
    bindingA = pd.read_csv("./deconv/data/Fig2A-FcgRI.csv")
    bindingB = pd.read_csv("./deconv/data/Fig2B-FcgRIIa-131H.csv")
    bindingC = pd.read_csv("./deconv/data/Fig2C-FcgRIIa-131R.csv")
    bindingD = pd.read_csv("./deconv/data/Fig2D-FcgRIIb.csv")
    bindingE = pd.read_csv("./deconv/data/Fig2E-FcgRIIIa-158F.csv")
    bindingF = pd.read_csv("./deconv/data/Fig2F-FcgRIIIa-158V.csv")
    bindingG = pd.read_csv("./deconv/data/Fig2G-FcgRIIIb-NA1.csv")
    bindingH = pd.read_csv("./deconv/data/Fig2H-FcgRIIIb-NA2.csv")

    mean_3a = (adcc3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc3b.groupby(level=0).sum()) / 4
    bindingAMean = np.nanmean(bindingA, axis=0)
    bindingBMean = np.nanmean(bindingB, axis=0)
    bindingCMean = np.nanmean(bindingC, axis=0)
    bindingDMean = np.nanmean(bindingD, axis=0)
    bindingEMean = np.nanmean(bindingE, axis=0)
    bindingFMean = np.nanmean(bindingF, axis=0)
    bindingGMean = np.nanmean(bindingG, axis=0)
    bindingHMean = np.nanmean(bindingH, axis=0)

    data = []
    data.append(mean_3a)
    data.append(mean_3b)
    data.append(bindingAMean)
    data.append(bindingBMean)
    data.append(bindingCMean)
    data.append(bindingDMean)
    data.append(bindingEMean)
    data.append(bindingFMean)
    data.append(bindingGMean)
    data.append(bindingHMean)

    data2 = np.transpose(np.array(data))

    pca2.fit(data2)
    data_new = pca2.transform(data2)

    ax[0].scatter(data_new[:,0], data_new[:,1])
    ax[0].set_title("Scores")
    ax[0].set_xlabel("Component 1 ({ratio:.2f})".format(ratio = pca2.explained_variance_ratio_[0]))
    ax[0].set_ylabel("Component 2 ({ratio:.2f})".format(ratio = pca2.explained_variance_ratio_[1]))

    l = ['adcc3a', 'adcc3b', 'bindingA', 'bindingB', 'bindingC', 'bindingD', 'bindingE', 'bindingF', 'bindingG', 'bindingH']
    loadings = pd.DataFrame(pca2.components_.T[:,:2], columns=['PC1', 'PC2'], index=l)
    ax[1].scatter(loadings.iloc[:,0], loadings.iloc[:,1])
    ax[1].set_title("Loadings")
    ax[1].set_xlabel("Component 1 ({ratio:.2f})".format(ratio = pca2.explained_variance_ratio_[0]))
    ax[1].set_ylabel("Component 2 ({ratio:.2f})".format(ratio = pca2.explained_variance_ratio_[1]))
    
    for i in range(10):
        ax[1].annotate(l[i], (loadings.iloc[i,0], loadings.iloc[i,1]))

    acc_variance = pca2.explained_variance_ratio_.copy()
    for i in range(1, 5):
        acc_variance[i] += acc_variance[i-1]
    ax[2].plot(range(1,6), acc_variance)
    ax[2].set_xticks(ticks = range(1,6))
    ax[2].set_ylabel("Explained Variance")
    ax[2].set_xlabel("Number of Components")

    # Add subplot labels
    subplotLabel(ax)

    return f