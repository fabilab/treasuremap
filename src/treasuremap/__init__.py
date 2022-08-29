import numpy as np
import igraph as ig

import _treasuremap


def layout_treasuremap(
        graph : ig.Graph,
        dist=None,
        seed=None,
        is_fixed=None,
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
    ):

    nvertices = graph.vcount()
    nedges = graph.ecount()

    if nvertices == 0:
        return ig.Layout([])
    if nvertices == 1:
        return ig.Layout([[0, 0]])

    if is_fixed is None:
        is_fixed = [False] * nvertices
    if len(is_fixed) != nvertices:
        raise ValueError("is_fixed must be a boolean vector of length n. vertices")

    edges = graph.get_edgelist()
    if dist is None:
        dist = [1] * nedges
    if len(dist) != nedges:
        raise ValueError("dist must be a vector with length n. edges")

    if min_dist < 0:
        raise ValueError(f"Minimum distance must be positive, got {min_dist}")

    if epochs < 0:
        raise ValueError(f"Number of epochs must be non-negative, got {min_dist}")

    if (sampling_prob <= 0) or (sampling_prob > 1):
        raise ValueError(f"Sampling probability must be in (0, 1], for {sampling_prob}")

    if dim not in (2, 3):
        raise ValueError(f"Number of dimensions must be 2 or 3, for {dim}")

    # Call C function
    result = _treasuremap.layout_treasuremap(
        nvertices, nedges,
        edges, dist, seed, is_fixed,
        min_dist, sampling_prob, epochs, dim,
    )

    #result = np.asarray(result).reshape((nvertices, 2))

    ## Recenter
    #result -= result.mean(axis=0)

    ## Make Layout
    #result = ig.Layout(list(result))
    return result
