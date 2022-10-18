# vim: fdm=indent
# author:     Fabio Zanini
# date:       5/08/19
# content:    Atlas subsampling
__all__ = ['subsample_atlas', 'average_atlas']

import warnings
import numpy as np
import pandas as pd
import anndata


def subsample_atlas(
        atlas,
        embedding_suffix='umap',
        n_grid_box_per_axis=30,
        ):
    '''Subsample atlas across cell types

    Args:
        atlas (anndata.Anndata): The dataset to subsample.
        embedding_suffix (str): Where to find the embedding metadata. If this
            argument is e.g. "umap", then treasuremap will look for the
            metadata in atlas.obsm["X_umap"].
        n_grid_box_per_axis (int): How fine grained the embedding grid should
            be. The default of 30 boxes means that a subsample based on a 2D
            embedding will contain at most 900 cells. Grid units with no cells
            will obviously be skipped.
    '''
    X_emb = atlas.obsm[f'X_{embedding_suffix}']
    ng = n_grid_box_per_axis  # Just s shorthand
    ndim = X_emb.shape[1]

    # Bin the data in each dimension
    cuts = []
    for i in range(ndim):
        cut_i = pd.cut(X_emb[:, i], bins=ng, labels=False, include_lowest=True)
        cuts.append(cut_i)
    cuts = pd.DataFrame(
        np.vstack(cuts).T, index=atlas.obs_names, columns=np.arange(ndim))
    cuts['_'] = 1

    # Group by box and pick cells
    cellnames = []
    for idx, rows in cuts.groupby(list(range(ndim))):
        # Pick at random within each grid box
        # NOTE: other algos are possible here, e.g. the closest to the box
        # center, or the closest to the group centroid
        cellname = rows.index[np.random.randint(len(rows))]
        # Append one cell to subsample
        cellnames.append(cellname)

    subsample = atlas[cellnames]

    return subsample
