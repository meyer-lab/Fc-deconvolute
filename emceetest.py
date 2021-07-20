#!/usr/bin/env python
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import std
from scipy.special import gamma, factorial
from deconv.imports import load_tables, load_figures, infer_x_fixed, ADCC_groups
from scipy.stats import norm
import emcee 

A_antiD, A_antiTNP, _, mixtures = load_tables()
adcc_3a, adcc_3b = load_figures()
mean_3a = (adcc_3a.groupby(level=0).sum()) / 4
mean_3b = (adcc_3b.groupby(level=0).sum()) / 4


def log_prob(x, mixtures, ADCC_true):
    k = 10
    x = np.transpose(x)
    ADCC_pred =  mixtures @ x
    residuals = ADCC_true - ADCC_pred
    scale = np.std(residuals)

    log_probability = norm.logpdf(residuals, loc=0, scale= scale)
    return log_probability

ndim, nwalkers = 24, 100
p0 = np.random.randn(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [A_antiD, mean_3a])
sampler.run_mcmc(p0, 10000)


# Traceback (most recent call last):
#   File "/home/ehunter/Fc-deconvolute/venv/lib/python3.9/site-packages/emcee/ensemble.py", line 619, in __call__
#     return self.f(x, *self.args, **self.kwargs)
#   File "/home/ehunter/Fc-deconvolute/emceetest.py", line 21, in log_prob
#     ADCC_pred =  x @ mixtures
# ValueError: matmul: Input operand 1 has a mismatch in its core dimension 0, with gufunc signature (n?,k),(k,m?)->(n?,m?) (size 20 is different from 24)
# Traceback (most recent call last):
# In[ ]:




