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
        # Make a DataFrame with cell type info and whether it's already picked
        picked = atlas.obs[[cell_type_column]].copy()
        picked['picked'] = False
        picked.loc[cellnames, 'picked'] = True
        picked['total'] = 1

        tmp = picked.groupby(cell_type_column).sum()
        cell_type_picked = tmp['picked']
        cell_type_total = tmp['total']
        cell_type_missing = minimum_cells_per_type - cell_type_picked

        gby = picked.groupby(['picked', cell_type_column])
        for cell_type in cell_type_picked.index:
            if cell_type_picked[cell_type] == cell_type_total[cell_type]:
                continue
            miss = cell_type_missing[cell_type]
            if miss <= 0:
                continue
            group = gby.get_group((False, cell_type)).index.values
            if miss < len(group):
                group = np.random.choice(group, size=miss, replace=False)
            picked.loc[group, 'picked'] = True

        cellnames = picked.loc[picked['picked']].index.values

    subsample = atlas[cellnames]

    return subsample
