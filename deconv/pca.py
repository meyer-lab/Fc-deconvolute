import numpy as np


def pca(X):
    """Fit the model by computing full SVD on X."""
    # Center data
    X -= np.mean(X, axis=0)

    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    # flip eigenvectors' sign to enforce deterministic output
    U, Vt = svd_flip(U, Vt)

    # Get variance explained by singular values
    explained_variance_ratio_ = np.square(S)
    explained_variance_ratio_ /= np.sum(explained_variance_ratio_)
    return U @ np.diag(S), Vt, explained_variance_ratio_


def svd_flip(u, v):
    """Sign correction to ensure deterministic output from SVD."""
    # columns of u, rows of v
    max_abs_cols = np.argmax(np.abs(u), axis=0)
    signs = np.sign(u[max_abs_cols, range(u.shape[1])])
    u *= signs
    v *= signs[:, np.newaxis]
    return u, v
