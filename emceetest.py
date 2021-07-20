import numpy as np
import matplotlib.pyplot as plt
from deconv.imports import load_tables, load_figures
from scipy.stats import norm
import emcee 

A_antiD, A_antiTNP, _, mixtures = load_tables()
adcc_3a, adcc_3b = load_figures()
mean_3a = (adcc_3a.groupby(level=0).sum()).to_numpy() / 4
mean_3b = (adcc_3b.groupby(level=0).sum()).to_numpy() / 4

def log_prob(x, mixtures, ADCC_true):
    residuals = ADCC_true - (mixtures @ x.T)
    logpdf = np.sum(norm.logpdf(residuals, loc=0, scale=1.0))

    # Enforce that inferred activities are positive
    bounds = 1000.0 * np.linalg.norm(np.maximum(x, 0.0))

    return logpdf - bounds

ndim, nwalkers = 24, 50
p0 = np.absolute(np.random.randn(nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [A_antiD, mean_3a])
sampler.run_mcmc(p0, 1000, progress=True, thin=5)
samples = sampler.get_chain(flat=True, discard=10)

plt.plot(samples[:, (9, 10, 11, 12, 13)])
plt.xlabel("step number")
plt.savefig("out.png")
