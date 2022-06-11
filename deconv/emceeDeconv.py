import numpyro
from numpyro.distributions import HalfNormal, Normal
from jax.random import PRNGKey
import jax.numpy as jnp
import arviz as az
from jax.config import config
from .imports import load_dekkers


config.update('jax_platforms', 'cpu')


def deconvModel(A_antiD, observed, error):
    activity = numpyro.sample("activity", HalfNormal(jnp.broadcast_to(5.0, (24, 12))))
    predict = jnp.dot(A_antiD, activity)

    # Use the mean of the sem across experimental replicates for error term
    sigm = jnp.broadcast_to(error, predict.shape)
    numpyro.sample("fit", Normal(predict, sigm), obs=observed)


def getEmceeTrace(return_az=False):
    data_dekkers = load_dekkers()

    A_antiD = data_dekkers["antiD"]
    X = data_dekkers["profiling"].values
    error = data_dekkers["profiling_std_mean"].values

    nuts_kernel = numpyro.infer.NUTS(deconvModel)
    mcmc = numpyro.infer.MCMC(nuts_kernel, num_warmup=500, num_samples=1000, num_chains=3, chain_method="vectorized")
    mcmc.run(PRNGKey(0), A_antiD, X, error)

    samples = mcmc.get_samples()
    data = az.from_numpyro(mcmc)
    assert jnp.amin(az.ess(data)["activity"].to_numpy()) > 100
    assert jnp.amax(az.rhat(data)["activity"].to_numpy()) < 1.01

    if return_az:
        return data

    return samples["activity"]
