import numpy as np

from _treasuremap import fit as _fit


def fit():

    edges = [[0, 1], [4, 5], [0, 3], [0, 2]]
    nvertices = 6

    res = np.zeros((nvertices, 2))
    edges = np.array(edges)
    distances = np.ones(len(edges))
    is_fixed = np.zeros(nvertices, int)
    min_dist = 0.01
    sampling_probability = 1.0
    epochs = 10
    ndim = 2

    return _fit(
        res,
        edges,
        distances,
        is_fixed,
        min_dist,
        sampling_probability,
        epochs,
        ndim,
    )
