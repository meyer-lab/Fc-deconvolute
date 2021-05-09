from sklearn.decomposition import PCA
from ..imports import load_tables
import matplotlib.pyplot as plt
from .common import subplotLabel, getSetup

ax, f = getSetup((6, 3), (1, 2))

A_antiD, A_antitnp, glycans, _ = load_tables()
pca = PCA(n_components=2)
pca.fit(A_antiD)
pca.fit(A_antitnp)

A_antiD_new = pca.transform(A_antiD)
A_antiD_new.shape
plt.scatter(A_antiD_new[:,0], A_antiD_new[:,1])
plt.title("PCA Plot of Glycan Matrix Scores")
plt.xlabel("Component 1")
plt.ylabel("Component 2")