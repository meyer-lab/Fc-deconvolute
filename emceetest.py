import numpy as np
import seaborn as sns
import corner
import matplotlib.pyplot as plt
from deconv.imports import load_tables, load_figures
from numba import njit
import emcee 

A_antiD, A_antiTNP, _, mixtures = load_tables()
adcc_3a, adcc_3b = load_figures()
mean_3a = (adcc_3a.groupby(level=0).sum()).to_numpy() / 4
mean_3b = (adcc_3b.groupby(level=0).sum()).to_numpy() / 4

@njit
def log_prob(x, mixtures, ADCC_true):
    residuals = ADCC_true - (mixtures @ x.T)
    residuals /= 0.5 # Scale residuals by stdev
    # Use expression for logpdf so we can compile with numba
    logpdf = np.sum(-np.square(residuals) / 2.0 - np.sqrt(2*np.pi))

    # Enforce that inferred activities are positive
    bounds = 1e6 * np.linalg.norm(np.minimum(x, 0.0))
    return logpdf - bounds

ndim, nwalkers = 24, 100
p0 = np.absolute(np.random.randn(nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [A_antiD, mean_3a])
sampler.run_mcmc(p0, 40000, progress=True)

burnin = 20000
thin = 10
samples = sampler.get_chain(discard=burnin, flat=True, thin=thin)
log_prob_samples = sampler.get_log_prob(discard=burnin, flat=True, thin=thin)

print("burn-in: {0}".format(burnin))
print("thin: {0}".format(thin))
print("flat chain shape: {0}".format(samples.shape))
print("flat log prob shape: {0}".format(log_prob_samples.shape))

all_samples = np.concatenate(
    (samples, log_prob_samples[:, None]), axis=1
)

labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
labels += ["log prob"]

sns.boxplot(data=samples, showfliers=False)
#corner.corner(all_samples, labels=labels);
plt.savefig('out.jpg')