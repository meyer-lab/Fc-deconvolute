from ..imports import load_tables, infer_x, load_bindingData, infer_x_fixed, ADCC_groups, R1_groups
from .common import subplotLabel, getSetup
import numpy as np

def makeFigure():
    ax, f = getSetup((6, 12), (4, 2))

    A_antiD, _, _, _ = load_tables()
    mean_binding = list(load_bindingData())
    R3_groups = ADCC_groups()
    R1_group = R1_groups()

    for ii, data in enumerate(mean_binding):
        mean_binding[ii] = mean_binding[ii].groupby(level=0).sum() / mean_binding[ii].groupby(level=0).count()

    infer_binding = []

    for x in range(1):
        for i in range(20):
            bind = mean_binding[x].drop(labels=i+1)
            A = np.delete(A_antiD, i, axis=0)

            glycans = infer_x_fixed(A, bind, R1_group)

            infer_binding.append(A_antiD[i] @ glycans)

        ax[x].scatter(mean_binding[x], infer_binding)
        ax[x].set_xlabel("Measured Binding")
        ax[x].set_ylabel("Inferred Binding")

        infer_binding = []

    for x in range(1,4):
        for i in range(20):
            bind = mean_binding[x].drop(labels=i+1)
            A = np.delete(A_antiD, i, axis=0)

            glycans = infer_x(A, bind)

            infer_binding.append(A_antiD[i] @ glycans)

        ax[x].scatter(mean_binding[x], infer_binding)
        ax[x].set_xlabel("Measured Binding")
        ax[x].set_ylabel("Inferred Binding")

        infer_binding = []

    for x in range(4,8):
        for i in range(20):
            bind = mean_binding[x].drop(labels=i+1)
            A = np.delete(A_antiD, i, axis=0)

            glycans = infer_x_fixed(A, bind, R3_groups)

            infer_binding.append(A_antiD[i] @ glycans)

        ax[x].scatter(mean_binding[x], infer_binding)
        ax[x].set_xlabel("Measured Binding")
        ax[x].set_ylabel("Inferred Binding")

        infer_binding = []

    ax[0].plot(range(3), range(3), color='r', linestyle='dashed')
    ax[1].plot(range(3), range(3), color='r', linestyle='dashed')
    ax[2].plot(range(3), range(3), color='r', linestyle='dashed')
    ax[3].plot(range(4), range(4), color='r', linestyle='dashed')
    ax[4].plot(range(70), range(70), color='r', linestyle='dashed')
    ax[5].plot(range(50), range(50), color='r', linestyle='dashed')
    ax[6].plot(range(30), range(30), color='r', linestyle='dashed')
    ax[7].plot(range(20), range(20), color='r', linestyle='dashed')

    ax[0].set_title("FcγRI")
    ax[1].set_title("FcγRIIa-131H")
    ax[2].set_title("FcγRIIa-131R")
    ax[3].set_title("FcγRIIb/c")
    ax[4].set_title("FcγRIIIa-158F")
    ax[5].set_title("FcγRIIIa-158V")
    ax[6].set_title("FcγRIIIb-NA1")
    ax[7].set_title("FcγRIIIb-NA2")

    subplotLabel(ax)

    return f
