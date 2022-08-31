# vim: fdm=indent
'''
author:     Fabio Zanini
date:       31/08/22
content:    Test algo on fly atlas data.
'''
import sys
import numpy as np
import pandas as pd
import anndata
import treasuremap
#import northstar

import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':


    print('Load data')
    fn_data = '/home/fabio/university/PI/projects/flyageingatlas/data/example/TS_Lung.h5ad'
    adata = anndata.read_h5ad(fn_data)

    print('Prepare partial embedding')
    nfixed = 0
    seed_name = 'umap_partial'
    adata.obsm[f'X_{seed_name}'] = adata.obsm['X_umap'].copy()
    adata.obsm[f'X_{seed_name}'][nfixed:] = 0
    is_fixed = np.ones(adata.shape[0], bool)
    is_fixed[nfixed:] = False

    print('Run treasuremap')
    treasuremap.treasuremap_adata(
        adata,
        seed_name=seed_name,
        is_fixed=is_fixed,
        epochs=100,
    )

    print('Plot embeddings')
    cell_types = adata.obs['cell_ontology_class'].cat.categories
    cmap = dict(zip(cell_types, sns.color_palette('husl', n_colors=len(cell_types))))
    gby = adata.obs.groupby('cell_ontology_class')
    fig, axs = plt.subplots(1, 2, figsize=(8, 4))
    for cell_type, group in gby:
        cellnames = group.index
        x, y = adata[cellnames].obsm['X_umap'].T
        axs[0].scatter(x, y, color=cmap[cell_type], alpha=0.1)

        x, y = adata[cellnames].obsm['X_treasuremap'].T
        if nfixed:
            x_fixed, y_fixed = x[:nfixed], y[:-nfixed]
            axs[1].scatter(x_fixed, y_fixed, color=cmap[cell_type], alpha=0.1)
        x_new, y_new = x[nfixed:], y[nfixed:]
        axs[1].scatter(x_new, y_new, marker='s', color=cmap[cell_type], alpha=0.1)
    axs[0].grid(True)
    axs[1].grid(True)
    fig.tight_layout()

    plt.ion(); plt.show()
