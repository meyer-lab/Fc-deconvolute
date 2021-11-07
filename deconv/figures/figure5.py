import numpy as np
from .common import subplotLabel, getSetup
from ..imports import infer_x, ADCC_groups, R1_groups, infer_x_fixed, load_dekkers
from sklearn.utils import resample


def makeFigure():
    ax, f = getSetup((6, 12), (4, 2))

    data_dekkers = load_dekkers()

    mean_binding = list(data_dekkers["bindings"])
    binding = data_dekkers["bindings"]
    A_antiD, glycan_list = data_dekkers["antiD"], data_dekkers["glycans"]
    R3_groups = ADCC_groups()
    R1_group = R1_groups()

    for ii, data in enumerate(mean_binding):
        mean_binding[ii] = mean_binding[ii].groupby(level=0).sum() / mean_binding[ii].groupby(level=0).count()

    deconvA_grouped = infer_x_fixed(A_antiD, mean_binding[0], R1_group)
    deconvB = infer_x(A_antiD, mean_binding[1])
    deconvC = infer_x(A_antiD, mean_binding[2])
    deconvD = infer_x(A_antiD, mean_binding[3])
    deconvE_grouped = infer_x_fixed(A_antiD, mean_binding[4], R3_groups)
    deconvF_grouped = infer_x_fixed(A_antiD, mean_binding[5], R3_groups)
    deconvG_grouped = infer_x_fixed(A_antiD, mean_binding[6], R3_groups)
    deconvH_grouped = infer_x_fixed(A_antiD, mean_binding[7], R3_groups)
    deconv = [deconvA_grouped, deconvB, deconvC, deconvD, deconvE_grouped, deconvF_grouped, deconvG_grouped, deconvH_grouped]

    glycans = []

    num_iters = 100

    for ii in range(1):
        for i in range(num_iters):
            new_bind = resample(binding[ii], n_samples=60, replace=True, stratify=binding[ii].index)
            mean = new_bind.groupby(level=0).sum() / new_bind.groupby(level=0).count()

            glycans.append(infer_x_fixed(A_antiD, mean, R1_group))

        new_glycans = np.array(glycans)
        error = (abs(deconv[ii] - np.quantile(new_glycans, 0.33, axis=0)), abs(np.quantile(new_glycans, 0.67, axis=0) - deconv[ii]))
        ax[ii].bar(glycan_list, deconv[ii], yerr=error)
        glycans = []

    for ii in range(1, 4):
        for i in range(num_iters):
            new_bind = resample(binding[ii], n_samples=60, replace=True, stratify=binding[ii].index)
            mean = new_bind.groupby(level=0).sum() / new_bind.groupby(level=0).count()

            glycans.append(infer_x(A_antiD, mean))

        new_glycans = np.array(glycans)
        error = (abs(deconv[ii] - np.quantile(new_glycans, 0.33, axis=0)), abs(np.quantile(new_glycans, 0.67, axis=0) - deconv[ii]))
        ax[ii].bar(glycan_list, deconv[ii], yerr=error)
        glycans = []

    for ii in range(4, 8):
        for i in range(num_iters):
            new_bind = resample(binding[ii], n_samples=60, replace=True, stratify=binding[ii].index)
            mean = new_bind.groupby(level=0).sum() / new_bind.groupby(level=0).count()

            glycans.append(infer_x_fixed(A_antiD, mean, R3_groups))

        new_glycans = np.array(glycans)
        error = (abs(deconv[ii] - np.quantile(new_glycans, 0.33, axis=0)), abs(np.quantile(new_glycans, 0.67, axis=0) - deconv[ii]))
        ax[ii].bar(glycan_list, deconv[ii], yerr=error)
        glycans = []

    ax[0].set_title("FcγRI")
    ax[1].set_title("FcγRIIa-131H")
    ax[2].set_title("FcγRIIa-131R")
    ax[3].set_title("FcγRIIb/c")
    ax[4].set_title("FcγRIIIa-158F")
    ax[5].set_title("FcγRIIIa-158V")
    ax[6].set_title("FcγRIIIb-NA1")
    ax[7].set_title("FcγRIIIb-NA2")

    for x in range(8):
        ax[x].set_xlabel("Glycans")
        ax[x].set_ylabel("Receptor Binding")
        ax[x].set_xticklabels(glycan_list, rotation=90)

    subplotLabel(ax)

    return f
