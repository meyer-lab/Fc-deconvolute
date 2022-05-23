import numpyro
from numpyro.distributions import HalfNormal, Normal
from jax.random import PRNGKey
import jax.numpy as jnp
from .imports import load_dekkers


def deconvModel(A_antiD, observed):
    activity = numpyro.sample('activity', HalfNormal(jnp.broadcast_to(5.0, (24, 12))))
    predict = jnp.dot(A_antiD, activity)
    numpyro.sample('fit', Normal(0.0, 0.1), obs=predict - observed)


def getEmceeTrace():
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]
    X = data_dekkers["profiling"].values

    nuts_kernel = numpyro.infer.NUTS(deconvModel)
    mcmc = numpyro.infer.MCMC(nuts_kernel, num_warmup=500, num_samples=1500, num_chains=4, chain_method="vectorized")
    mcmc.run(PRNGKey(0), A_antiD, X)

    samples = mcmc.get_samples()
    summary = numpyro.diagnostics.summary(mcmc.get_samples(group_by_chain=True))["activity"]
    assert jnp.amax(summary["r_hat"]) < 1.01
    assert jnp.amax(summary["n_eff"]) > 1000
    return samples["activity"]
