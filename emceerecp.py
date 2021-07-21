import numpy as np
from scipy.stats._continuous_distns import _norm_logpdf
import seaborn as sns
import corner
import matplotlib.pyplot as plt
from deconv.imports import load_tables, load_bindingData
import emcee

A_antiD, A_antiTNP, _, mixtures = load_tables()

mean_binding = list(load_bindingData())
mean_binding = np.array([m.groupby(level=0).mean() for m in mean_binding])


def log_prob(x, mixtures, trueVals):
    x_x = np.reshape(x[:, 0:72], (x.shape[0], 24, 3))
    x_y = np.reshape(x[:, 72::], (x.shape[0], 3, 8))
    np.testing.assert_array_equal(x[:, 0], x_x[:, 0, 0])
    XX = np.einsum('ijk,ikl->ijl', x_x, x_y)

    residuals = trueVals[:, :, np.newaxis] - (mixtures @ XX).T
    logpdf = np.sum(_norm_logpdf(residuals / 1.0), axis=(0, 1))

    # Enforce scale
    x_scale = np.linalg.norm(x_x, axis=1)
    scale_pdf = np.sum(_norm_logpdf((x_scale - 1.0) * 100.0), axis=1)

    # Enforce that inferred activities are positive
    bounds = 1e4 * np.linalg.norm(np.minimum(x, 0.0), axis=1)
    return logpdf - bounds + scale_pdf


ndim, nwalkers = 96, 200
p0 = np.absolute(np.random.randn(nwalkers, ndim))

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, vectorize=True, args=[A_antiD, mean_binding])
sampler.run_mcmc(p0, 50000, progress=True)

burnin = 20000
thin = 20
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

sns.boxplot(data=samples[:, 0:72], showfliers=False)
#plt.plot(log_prob_samples)
#corner.corner(all_samples, labels=labels);
plt.savefig('out.jpg')
