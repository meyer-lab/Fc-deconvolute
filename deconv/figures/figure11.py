import numpy as np
import matplotlib.pyplot as plt
from .common import getSetup
from ..emceeDeconv import getEmceeTrace

def makeFigure():
    ax, f = getSetup((12, 6), (1, 3))

    trace = getEmceeTrace()

    activity_scores = trace.posterior.activity_scores[0]
    median = np.median(activity_scores, axis=0)
    p33 = np.percentile(activity_scores, 33, axis=0)
    p66 = np.percentile(activity_scores, 66, axis=0)

    plt.rcParams["figure.figsize"] = (20,3)

    for i in range(3):
        ax[i].set_title("Activity Scores [{}]".format(i))
        ax[i].errorbar(range(1,25), median[:,i], yerr=[p33[:,i], p66[:,i]], fmt='o')

    return f
