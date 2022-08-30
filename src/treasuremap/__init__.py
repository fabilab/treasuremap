try:
    from igraph import Graph, Layout
except ImportError:
    Graph = None
    Layout = None
try:
    from anndata import AnnData
except ImportError:
    AnnData = None

import _treasuremap


def _recenter(result):
    """Recenter layout after algo"""
    dim = len(result[0])

    # Find center
    centers = [0 for i in range(dim)]
    for cell in result:
        for i, coord in enumerate(cell):
            centers[i] += coord
    for i in range(dim):
        centers[i] /= len(result)

    # Recenter in place
    for cell in result:
        for i, coord in enumerate(cell):
            cell[i] -= centers[i]


def _layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist=None,
        seed=None,
        is_fixed=None,
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
    ):

    if min_dist < 0:
        raise ValueError(f"Minimum distance must be positive, got {min_dist}")

    if epochs < 0:
        raise ValueError(f"Number of epochs must be non-negative, got {min_dist}")

    if (sampling_prob <= 0) or (sampling_prob > 1):
        raise ValueError(f"Sampling probability must be in (0, 1], for {sampling_prob}")

    if dim not in (2, 3):
        raise ValueError(f"Number of dimensions must be 2 or 3, for {dim}")

    if nvertices == 0:
        return Layout([])
    if nvertices == 1:
        return Layout([[0, 0]])

    if len(edges) == 0:
        raise ValueError("graph has no edges, Treasuremap cannot do much...")

    if is_fixed is None:
        is_fixed = [False] * nvertices
    if len(is_fixed) != nvertices:
        raise ValueError("is_fixed must be a boolean vector of length n. vertices")

    if dist is None:
        dist = [1] * nedges
    if len(dist) != nedges:
        raise ValueError("dist must be a vector with length n. edges")

    # Call C function
    result = _treasuremap.layout_treasuremap(
        nvertices, nedges,
        edges, dist, seed, is_fixed,
        min_dist, sampling_prob, epochs, dim,
    )

    # Recenter
    _recenter(result)

    # Make Layout
    result = Layout(list(result))
    return result


def treasuremap_adata(
        adata : AnnData,
        seed_name='UMAP',
        is_fixed=None,
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
    ):

    nvertices = adata.shape[0]

    dist_matrix = adata.obsp['distances'].tocoo()
    edges = list(zip(dist_matrix.rows, dist_matrix.cols))
    dist = dist_matrix.data
    nedges = len(dist)

    seed = adata.obsm[f'X_{seed_name}'].tolist()

    return _layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist,
        seed,
        is_fixed,
        min_dist,
        sampling_prob,
        epochs,
        dim,
    )


def treasuremap_igraph(
        graph : Graph,
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
    edges = graph.get_edgelist()

    result = _layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist,
        seed,
        is_fixed,
        min_dist,
        sampling_prob,
        epochs,
        dim,
    )

    result = Layout(result)
    return result
