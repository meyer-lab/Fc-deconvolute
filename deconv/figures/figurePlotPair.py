import arviz as az
import matplotlib.pyplot as plt
from ..emceeDeconv import getEmceeTrace


def makeFigure():
    activity = getEmceeTrace(return_az=True)

    az.plot_pair(activity)
    plt.savefig('output/cornerplot.png')

    return None
