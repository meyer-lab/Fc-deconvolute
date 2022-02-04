import numpy as np
from .common import getSetup, subplotLabel
from ..emceeDeconv import getEmceeTrace
from ..imports import load_dekkers

def makeFigure():
    ax, f = getSetup((7.5,7.5), (4, 3))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]
    l = ['Binding FcγRIa', 'Binding FcγRIIa 131H', 'Binding FcγRIIa 131R', 'Binding FcγRIIb/c', 'Binding FcγRIIIa 158F', 'Binding FcγRIIIa 158V', 'Binding Fc-FcγRIIIb NA1', 'Binding Fc-FcγRIIIb NA2', 'ADCC FcγRIIIA158F/F', 'ADCC FcγRIIIA158V/V', 'Complement Activation C1q', 'Complement Activation C4']

    trace = getEmceeTrace()

    activity_scores = trace.posterior.activity_scores[0]
    activity_loadings = trace.posterior.activity_loadings[0]

    full = np.einsum("ijk,ikl->ijl", activity_scores, activity_loadings)

    qqs = np.quantile(full, (0.33, 0.5, 0.66), axis=0)
    median = qqs[1, :, :]
    p33 = median - qqs[0, :, :]
    p66 = qqs[2, :, :] - median

    print(median.shape)
    print(len(glycans))

    for i in range(median.shape[1]):
        ax[i].errorbar(glycans, median[:, i], yerr=[p33[:, i], p66[:, i]], fmt='o', markersize = 4)
        ax[i].set_ylabel(l[i], size=6)
        ax[i].set_xticklabels(glycans, rotation=90, size=6)
        _,_,_,y2 = ax[i].axis()
        ax[i].set_ylim([0, y2]) 

    subplotLabel(ax)

    return f
