try:
    from igraph import Graph, Layout
except ImportError:
    Graph = None
    Layout = None
try:
    from anndata import AnnData
    import numpy as np
except ImportError:
    AnnData = None
    np = None

import treasuremap._treasuremap as _treasuremap


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
        seed_name='umap',
        is_fixed=None,
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
        copy=False,
    ):
    if AnnData is None:
        raise ImportError("Install the package anndata to use this function")

    nvertices = adata.shape[0]

    if 'distances' not in adata.obsp:
        raise KeyError("AnnData object must have an obsp['distances'] matrix")
    dist_matrix = adata.obsp['distances'].tocoo()

    edges = list(zip(dist_matrix.row, dist_matrix.col))
    dist = dist_matrix.data
    nedges = len(dist)

    obsm_name = f'X_{seed_name}'
    if obsm_name not in adata.obsm:
        raise KeyError("AnnData object missing {obsm_name} field")
    seed = adata.obsm[f'X_{seed_name}']
    seed_dim = seed.shape[1]
    if seed_dim != dim:
        raise ValueError("{obsm_name} is {seed_dim}D, requested {dim}D embedding")
    seed = seed.tolist()

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
    result = np.asarray(result).astype(np.float32)

    if copy:
        return result

    adata.obsm['X_treasuremap'] = result


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

    if Graph is None:
        raise ImportError("Install the package igraph to use this function")

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
