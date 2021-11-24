import numpy as np
import matplotlib.pyplot as plt
from .common import getSetup
from ..emceeDeconv import getEmceeTrace
from ..imports import load_dekkers

def makeFigure():
    ax, f = getSetup((12,12), (4, 3))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]
    l = ['Binding FcγRIa', 'Binding FcγRIIa 131H', 'Binding FcγRIIa 131R', 'Binding FcγRIIb/c', 'Binding FcγRIIIa 158F', 'Binding FcγRIIIa 158V', 'Binding Fc-FcγRIIIb NA1', 'Binding Fc-FcγRIIIb NA2', 'ADCC FcγRIIIA158F/F', 'ADCC FcγRIIIA158V/V', 'Complement Activation C1q', 'Complement Activation C4']

    
    trace = getEmceeTrace()

    activity_scores = trace.posterior.activity_scores[0]
    activity_loadings = trace.posterior.activity_loadings[0]

    median_scores = np.median(activity_scores, axis=0) 
    p33_scores = np.percentile(activity_scores, 33, axis=0)
    p66_scores = np.percentile(activity_scores, 66, axis=0)

    median_loadings = np.median(activity_loadings, axis=0)
    p33_loadings = np.percentile(activity_loadings, 33, axis=0)
    p66_loadings = np.percentile(activity_loadings, 66, axis=0)

    measurement_glycan_matrix = np.transpose(np.dot(median_scores, median_loadings))
    p33 = np.transpose(np.dot(p33_scores, p33_loadings))
    p66 = np.transpose(np.dot(p66_scores, p66_loadings))


    for i in range(len(measurement_glycan_matrix)):
        ax[i].errorbar(glycans, measurement_glycan_matrix[i], yerr=[p33[i], p66[i]], fmt='o')
        ax[i].set_title(l[i])
        ax[i].set_xticklabels(glycans, rotation=90)
