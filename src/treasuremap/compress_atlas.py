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
        minimum_cells_per_type=0,
        cell_type_column=None,
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
        minimum_cells_per_type (int): At the end of the subsampling, add back
            at least this number of cells for each annotated type. This option
            is useful if you want to make sure rare populations are represented
            even though they are contained in a small area of the embedding.
        cell_type_column (str or None): This option is only used if
            minimum_cells_per_type > 0. If that's the case, this is the name
            of the column in atlas.obs containing the cell type annotations.
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

    if minimum_cells_per_type > 0:
        cn_set = set(cellnames)
        cell_type_n = atlas.obs.loc[cellnames, cell_type_column].value_counts()
        cell_types_incomplete = cell_type_n.index[cell_type_n < minimum_cells_per_type]
        gby = atlas.obs.groupby(cell_type_column)
        for cell_type in cell_types_incomplete:
            nmiss = minimum_cells_per_type - cell_type_n[cell_type]
            cellnames_type = gby.get_group(cell_type).index
            cellnames_free = list(set(cellnames_type) - cn_set)
            if nmiss > len(cellnames_free):
                cellnames_take = cellnames_free
            else:
                cellnames_take = np.random.choice(
                        cellnames_free, size=nmiss, replace=True)
            cellnames.extend(cellnames_take)

    subsample = atlas[cellnames]

    return subsample
