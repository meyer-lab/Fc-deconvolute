from .common import getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers
import seaborn as sns

def makeFigure():
    ax, f = getSetup((7.5, 15), (6, 2))

    data_dekkers = load_dekkers()
    glycans = data_dekkers["glycans"]

    data = data_dekkers["profiling"]
    data2 = data.pivot(index="index", columns="receptor", values="binding")

    trace = getEmceeTrace()

    activity = trace.posterior.activity[0]

    qqs = np.quantile(activity, (0.025, 0.5, 0.975), axis=0)
    median = qqs[1, :, :]

    l = data2.columns

    d = {
    "Fucosylation": pd.Series(['Fucosylated','Fucosylated','Fucosylated',
        'Fucosylated','Fucosylated','Fucosylated','Fucosylated','Fucosylated',
        'Fucosylated','Fucosylated','Fucosylated','Fucosylated',"Afucosylated",
        "Afucosylated","Afucosylated","Afucosylated","Afucosylated","Afucosylated",
        "Afucosylated","Afucosylated","Afucosylated","Afucosylated","Afucosylated",
        "Afucosylated",]),
    "Galactosylation": pd.Series(["No Galactose","1 Galactose","2 Galactose","No Galactose",
        "1 Galactose","2 Galactose","1 Galactose","2 Galactose","2 Galactose","1 Galactose",
        "2 Galactose","2 Galactose","No Galactose","1 Galactose","2 Galactose","No Galactose",
        "1 Galactose","2 Galactose","1 Galactose","2 Galactose","2 Galactose","1 Galactose",
        "2 Galactose","2 Galactose"]),
    "Bisection": pd.Series(["None", "None", "None", "Bisection", "Bisection", "Bisection", 
        "None", "None", "None", "Bisection", "Bisection", "Bisection", "None", "None", "None", 
         "Bisection", "Bisection", "Bisection", "None", "None", "None","Bisection", "Bisection",
          "Bisection"]),
    "Sialylation": pd.Series(["No Sialic Acid", "No Sialic Acid", "No Sialic Acid", "No Sialic Acid",
         "No Sialic Acid", "No Sialic Acid", "1 Sialic Acid","1 Sialic Acid", "2 Sialic Acid", 
         "1 Sialic Acid","1 Sialic Acid", "2 Sialic Acid", "No Sialic Acid", "No Sialic Acid",
          "No Sialic Acid", "No Sialic Acid", "No Sialic Acid", "No Sialic Acid", "1 Sialic Acid",
          "1 Sialic Acid", "2 Sialic Acid", "1 Sialic Acid","1 Sialic Acid", "2 Sialic Acid" ]),
    "Sial Markers": pd.Series(['o','o','o','o','o','o','s', 's', 'v','s', 's', 'v','o','o','o','o','o','o', 's', 's', 'v','s', 's', 'v'])}
    df = pd.DataFrame(d)
    df["glycans"] = glycans

    for i in range(12):
        header = '{}'.format(l[i])
        df[header] = median[:,i]

    sns.set_theme(style="whitegrid", palette= ['dodgerblue','salmon'])
    for i in range(12):
        ax[i] = sns.scatterplot(ax = ax[i], data = df, x='Sialylation', y=l[i], hue='Fucosylation', 
            style = 'Bisection', markers=['d','X'])

    return(f)
