from .common import subplotLabel, getSetup
from ..emceeDeconv import getEmceeTrace
import pandas as pd
import numpy as np
from ..imports import load_dekkers
import seaborn as sns

def makeFigure():
    ax, f = getSetup((15,10), (4,3))

    #imports and formats
    trace = getEmceeTrace()
    data_dekkers = load_dekkers()
    l = ['Binding FcγRIa', 'Binding FcγRIIa 131H', 'Binding FcγRIIa 131R', 'Binding FcγRIIb/c', 'Binding FcγRIIIa 158F', 'Binding FcγRIIIa 158V', 'Binding Fc-FcγRIIIb NA1', 'Binding Fc-FcγRIIIb NA2', 'ADCC FcγRIIIA158F/F', 'ADCC FcγRIIIA158V/V', 'Complement Activation C1q', 'Complement Activation C4']

    activity_scores = trace.posterior.activity_scores[0]
    activity_loadings = trace.posterior.activity_loadings[0]
    full = np.einsum("ijk,ikl->ijl", activity_scores, activity_loadings)

    qqs = np.quantile(full, (0.33, 0.5, 0.66), axis=0)
    median = qqs[1, :, :]
    glycans = data_dekkers["glycans"]

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
          "1 Sialic Acid", "2 Sialic Acid", "1 Sialic Acid","1 Sialic Acid", "2 Sialic Acid" ])}
    df = pd.DataFrame(d)
    df["glycans"] = glycans

    for i in range(12):
        header = '{}'.format(l[i])
        df[header] = median[:,i]

    sns.set_theme(style="whitegrid", palette="muted")
    for i in range(12):
        ax[i] = sns.swarmplot(ax = ax[i], data=df, x="Galactosylation", y= l[i], hue = "Fucosylation")

    return(f)
