import numpy as np
import igraph as ig

from _treasuremap import _layout_treasuremap


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

    use_seed = seed is not None
    nvertices = graph.vcount()
    edges = graph.get_edgelist()
    edges = np.asarray(edges)

    if is_fixed is None:
        is_fixed = np.zeros(nvertices, int)
    elif isinstance(is_fixed, str):
        is_fixed = graph.vs[is_fixed]
    is_fixed = np.asarray(is_fixed)

    if dist is None:
        dist = np.ones(len(edges))

    if min_dist < 0:
        raise ValueError(f"Minimum distance must be positive, got {min_dist}")

    if epochs < 0:
        raise ValueError(f"Number of epochs must be non-negative, got {min_dist}")

    if (sampling_prob <= 0) or (sampling_prob > 1):
        raise ValueError(f"Sampling probability must be in (0, 1], for {sampling_prob}")

    if dim not in (2, 3):
        raise ValueError(f"Number of dimensions must be 2 or 3, for {dim}")

    # Call C function
    result = _layout_treasuremap(
        edges,
        dist,
        seed,
        use_seed,
        is_fixed,
        min_dist,
        sampling_prob,
        epochs,
        dim,
    )
    result = result.reshape((nvertices, 2))

    # Recenter
    result -= result.mean(axis=0)

    # Make Layout
    result = ig.Layout(list(result))
    return result
