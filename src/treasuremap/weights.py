# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/11/22
content:    Compute UMAP-like weights from distances.
'''
from treasuremap._treasuremap import compute_weights as _compute_weights


def compute_weights(
    distance_matrix,
    remove_zero_weight_edges=True,
    ):
    '''Compute UMAP-like weights from distance matrix

    Args:
        distance_matrix: sparse matrix containing distances associated with each
          edge. This matrix needs not be symmetric.
        remove_zero_weight_edges: whether to remove zero-weight edges from the
          returned matrix, which is symmetrized and triangular. Default is True.

    Returns: CSR sparse matrix with weights. It is triangular to save space, with
        zero-diagonal elements and an understanding that the weights are to be
        interpreted as symmetric. You can plug in this matrix directly into
        treasuremap for layouting.
    '''
    distance_matrix = distance_matrix.tocoo()

    nvertices = distance_matrix.shape[0]
    edges = list(zip(distance_matrix.row, distance_matrix.col))
    nedges = len(edges)
    distances = list(distance_matrix.data)

    weights = _compute_weights(
        nvertices, nedges, edges, distances,
    )

    weight_matrix = distance_matrix.copy()
    for i, weight in enumerate(weights):
        weight_matrix.data[i] = weight

    if remove_zero_weight_edges:
        weight_matrix.eliminate_zeros()

    weight_matrix = weight_matrix.tocsr()

    return weight_matrix
