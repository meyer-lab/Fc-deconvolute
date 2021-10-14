from deconv.imports import load_tables, load_figures, infer_x
from deconv.figures.common import subplotLabel, getSetup
import numpy as np

def makeFigure():
    ax, f = getSetup((6, 3), (1, 2))

    A_antiD, _, _, _ = load_tables()
    adcc3a, adcc3b = load_figures()

    mean_3a = (adcc3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc3b.groupby(level=0).sum()) / 4

    infer_adcc_3a = []
    infer_adcc_3b = []

    for i in range(20):
        adccA = mean_3a.drop(labels=i+1)
        adccB = mean_3b.drop(labels=i+1)
        A = np.delete(A_antiD, i, axis=0)

        glycans_3a = infer_x(A, adccA)
        glycans_3b = infer_x(A, adccB)

        infer_adcc_3a.append(A_antiD[i] @ glycans_3a)
        infer_adcc_3b.append(A_antiD[i] @ glycans_3b)


    ind1 = [0, 5, 10, 15, 20, 25, 30, 35]
    ind2 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]

    ax[0].scatter(mean_3a, infer_adcc_3a)
    ax[0].set_title("Original and Inferred ADCC (Fig. 3A)")
    ax[0].set_xlabel("Actual ADCC")
    ax[0].set_ylabel("Inferred ADCC")
    ax[0].set_xticks(ind1)
    ax[0].set_yticks(ind1)
    ax[0].plot(range(35), range(35), color='r', linestyle='dashed')


    ax[1].scatter(mean_3b, infer_adcc_3b)
    ax[1].set_title("Original and Inferred ADCC (Fig. 3B)")
    ax[1].set_xlabel("Actual ADCC")
    ax[1].set_ylabel("Inferred ADCC")
    ax[1].set_xticks(ind2)
    ax[1].set_yticks(ind2)
    ax[1].plot(range(45), range(45), color='r', linestyle='dashed')

    return f