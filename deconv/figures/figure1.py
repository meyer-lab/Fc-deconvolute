from matplotlib import pyplot as plt
from common import load_tables, load_figures, infer_x

def makeFig1():
    A_antiD, A_antiTNP, glycan_list = load_tables()
    adcc_3a, adcc_3b = load_figures()

    mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
    mean_3b = (adcc_3b.groupby(level=0).sum()) / 4

    glycans_3a = infer_x(A_antiD, mean_3a)
    glycans_3b = infer_x(A_antiD, mean_3b)

    plot_3a = plt.figure()
    ax = plot_3a.add_axes([0,0,1,1])
    ax.bar(glycan_list, glycans_3a)
    plt.title('ADCC (Fig. 3A)')
    plt.xlabel('Glycans')
    plt.xticks(rotation=90)

    plot_3b = plt.figure()
    ax = plot_3b.add_axes([0,0,1,1])
    ax.bar(glycan_list, glycans_3b)
    plt.title('ADCC (Fig. 3B)')
    plt.xlabel('Glycans')
    plt.xticks(rotation=90)

    return (plot_3a, plot_3b)

