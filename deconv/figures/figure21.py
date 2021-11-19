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
    median_loadings = np.median(activity_loadings, axis=0)
    measurement_glycan_matrix = np.transpose(np.dot(median_scores, median_loadings))

    for i in range(len(measurement_glycan_matrix)):
        ax[i].bar(glycans, measurement_glycan_matrix[i])
        ax[i].set_title(l[i])
        ax[i].set_xticklabels(glycans, rotation=90)
