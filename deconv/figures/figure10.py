from matplotlib.pyplot import get
import numpy as np
import arviz as az
from ..emceeDeconv import getEmceeTrace

def makeFigure():
    trace = getEmceeTrace()
    f = az.plot_forest(trace)
    #plt.savefig('out1.jpg')
    return f