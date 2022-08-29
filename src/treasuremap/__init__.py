import numpy as np

from _treasuremap import _layout_treasuremap


def fit(
    graph=None,
    seed=None,
    is_fixed=None,
    distances=None,
    min_dist=0.01,
    sampling_probability=1.0,
    epochs=10,
    ndim=2,
):

    # FIXME
    if graph is None:
        import igraph as ig
        nvertices = 6
        graph = ig.Graph(
            nvertices,
            [[0, 1], [4, 5], [0, 3], [0, 2]],
        )
        seed = np.arange(nvertices * 2).reshape((nvertices, 2))
        is_fixed = [0, 1, 0, 1, 0, 1]

    nvertices = graph.vcount()
    edges = graph.get_edgelist()

    if is_fixed is None:
        is_fixed = np.zeros(nvertices, int)
    elif isinstance(is_fixed, str):
        is_fixed = graph.vs[is_fixed]

    edges = np.asarray(edges)
    is_fixed = np.asarray(is_fixed)

    if distances is None:
        distances = np.ones(len(edges))

    # Call C function
    result = _layout_treasuremap(
        seed,
        edges,
        distances,
        is_fixed,
        min_dist,
        sampling_probability,
        epochs,
        ndim,
    )

    result = result.reshape((nvertices, 2))
    return result
