# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/11/22
content:    Build k-nearest neighbor graph for treasuremap.
'''
try:
    import numpy as np
except ImportError:
    np = None

try:
    import anndata
except ImportError:
    anndata = None

try:
    import scanpy as sc
except ImportError:
    sc = None

try:
    import hnswlib
except:
    hnswlib = None

def build_knn(
    adata,
    adata_reference,
    space='pca',
    features=1000,
    n_pcs=50,
    n_internal_neighbors=14,
    n_external_neighbors=5,
    verbose=0,
    metric='euclidean',
    approximate=False,
    ):
    """Build k-nearest neighbor graph between a data set and a reference"""
    from scipy.sparse import lil_matrix

    if (np is None) or (anndata is None) or (sc is None) or (hnswlib is None):
        raise ImportError(
                "This function requires numpy, anndata, scanpy and hnswlib.")

    n = adata.shape[0]
    n_reference = adata_reference.shape[0]
    ntot = n_reference + n

    if verbose:
        print('Concatenate AnnData objects')
    adatam = anndata.concat(
            [adata_reference, adata],
            join='outer',
            )

    if adatam.X.max() > 50:
        if verbose:
            print('Take log1p')
        sc.pp.log1p(adatam)

    if space == 'pca':
        if features != 0:
            if verbose:
                print('Select features')
            if isinstance(features, int):
                sc.pp.highly_variable_genes(adatam, n_top_genes=features)
            else:
                adatam.var['highly_variable'] = False
                adatam.var.loc[features, 'highly_variable'] = True

            adatam.raw = adatam
            adatam = adatam[:, adatam.var.highly_variable]

        #if verbose:
        #    print('Scale gene expression to 10')
        #sc.pp.scale(adata, max_value=10)

        if verbose:
            print('Principal Component Analysis')
        sc.tl.pca(adatam, n_comps=n_pcs)

    if verbose:
        print('Build knn')

    if approximate:
        if verbose:
            print('Internal neighbors (approximate)')
        ef, M = 200, 48  # FIXME: check this
        metric_hnsw = 'cosine' if metric == 'correlation' else 'l2'
        ref_index = hnswlib.Index(space=metric_hnsw, dim=n_pcs)
        ref_index.init_index(
            max_elements=n,
            ef_construction=ef,
            M=M,
        )
        ref_index.add_items(adatam.obsm['X_pca'][n_reference:], np.arange(n) + n_reference)
        ref_index.set_ef(ef)  # FIXME: Probably redundant?
        idx_internal, dist_internal = ref_index.knn_query(
                adatam.obsm[f'X_{space}'][n_reference:], k=n_internal_neighbors + 1)
        # Exclude loops/identities
        idx_internal = idx_internal[:, 1:]
        dist_internal = dist_internal[:, 1:]
        # Sometimes the feature selection results in zero expression
        dist_internal = np.maximum(dist_internal, 0).astype(np.float32)
        idx_internal = idx_internal.astype(int)

    ## Incoming edges from reference
    #idx_incoming, dist_incoming = ref_index.knn_query(
    #        adatam.obsm[f'X_{space}'][:n_reference], k=n_external_neighbors)
    #dist_incoming = np.maximum(dist_incoming, 0)

    else:
        if verbose:
            print('Internal neighbors (exact)')
        from scipy.spatial.distance import cdist
        idx_internal = np.zeros((n, n_internal_neighbors), int)
        dist_internal = np.zeros((n, n_internal_neighbors), np.float32)
        x = adatam.obsm[f'X_{space}'][n_reference:]
        i = 0
        blk = 100
        while i < n:
            print(f'knn vertex: {i+1} / {ntot}', end='\r')
            y = x[i : i + blk]
            d = cdist(y, x, metric=metric)

            # FIXME: Exclude internal neighbors for testing
            #d[:, n_reference:] = 1000

            idxi = d.argsort(axis=1)[:, 1:n_internal_neighbors + 1]
            disti = np.take_along_axis(d, idxi, 1)
            # Offset them since their indices do not actually start from zero
            idx_internal[i : i+blk] = idxi + n_reference
            dist_internal[i: i+blk] = disti
            i += blk
        dist_internal = np.maximum(dist_internal, 0)

        tmp = idx_internal.argsort(axis=1)
        idx_internal = np.take_along_axis(idx_internal, tmp, 1)
        dist_internal = np.take_along_axis(dist_internal, tmp, 1)

    if approximate:
        if verbose:
            print('External neighbors (approximate)')
        ef, M = 200, 48
        metric_hnsw = 'cosine' if metric == 'correlation' else 'l2'
        ref_index = hnswlib.Index(space=metric_hnsw, dim=n_pcs)
        ref_index.init_index(
            max_elements=n_reference,
            ef_construction=ef,
            M=M,
        )
        ref_index.add_items(adata.obsm[f'X_{space}'][:n_reference], np.arange(n_reference))
        ref_index.set_ef(ef)
        idx_external, dist_external = ref_index.knn_query(
                adatam.obsm[f'X_{space}'][n_reference:], k=n_external_neighbors)
        dist_external = np.maximum(dist_external, 0).astype(np.float32)
        idx_external = idx_external.astype(int)

    else:
        if verbose:
            print('External neighbors (exact)')
        from scipy.spatial.distance import cdist
        idx_external = np.zeros((n, n_external_neighbors), int)
        dist_external = np.zeros((n, n_external_neighbors), np.float32)
        x1 = adatam.obsm[f'X_{space}'][n_reference:]
        x2 = adatam.obsm[f'X_{space}'][:n_reference]
        i = 0
        blk = 100
        while i < n:
            print(f'knn vertex: {i+1} / {n}', end='\r')
            y = x1[i : i + blk]
            d = cdist(y, x2, metric=metric)

            idxi = d.argpartition(n_external_neighbors)[:, :n_external_neighbors]
            disti = np.take_along_axis(d, idxi, 1)
            idx_external[i: i + len(idxi)] = idxi
            dist_external[i: i + len(idxi)] = disti
            i += blk
        dist_external = np.maximum(dist_external, 0)

        tmp = idx_external.argsort(axis=1)
        idx_external = np.take_along_axis(idx_external, tmp, 1)
        dist_external = np.take_along_axis(dist_external, tmp, 1)

    if verbose:
        print('Merge external and internal neighbors')
        n_neis = n_internal_neighbors + n_external_neighbors
        idx = np.vstack([
            np.zeros((n_reference, n_neis), int),
            np.hstack([idx_internal, idx_external]),
            ])
        dist = np.vstack([
            np.zeros((n_reference, n_neis), np.float32),
            np.hstack([dist_internal, dist_external]),
            ])

        # Sort them for convenience
        tmp = dist.argsort(axis=1)
        idx = np.take_along_axis(idx, tmp, 1)
        dist = np.take_along_axis(dist, tmp, 1)

    ## FIXME
    #return {
    #    'obs_names': adatam.obs_names,
    #    'obs_names_free': adatam.obs_names[n_reference:],
    #    'cellnames_neighbors': adatam.obs_names[idx_external],
    #    'dist_external': dist_external,
    #}

    ##if verbose:
    ##    print('Perform connectivity correction')
    ##dist_int_min = np.percentile(dist_internal, 1)
    ##dist_ext_min = np.percentile(dist_external, 1)
    ##dist_in_min = np.percentile(dist_incoming, 1)
    ##if verbose:
    ##    print('dist_int_min: {:.3f}'.format(dist_int_min))
    ##    print('dist_ext_min: {:.3f}'.format(dist_ext_min))
    ##    print('dist_in_min: {:.3f}'.format(dist_in_min))
    ##if dist_ext_min > dist_int_min:
    ##    dist_external += dist_int_min - dist_ext_min
    ##    dist_external = np.maximum(dist_external, dist_int_min)
    ##if dist_in_min > dist_int_min:
    ##    dist_incoming += dist_int_min - dist_in_min
    ##    dist_incoming = np.maximum(dist_incoming, dist_int_min)

    ##if verbose:
    ##    print('Combine internal and external neighbors')
    ##dist = np.hstack([dist_internal, dist_external])
    ##idx = np.hstack([idx_internal, idx_external])
    ### Partition the closest neighbors from the rest
    ### NOTE: the neighbors themselves are unsorted, that should be fine
    ### since the function that computes rho/sigma does NOT assume the
    ### first neighbor is the closest one
    ###part = dist.argpartition(n_internal_neighbors)[:, :n_internal_neighbors]
    ###dist = np.take_along_axis(dist, part, 1)
    ###idx = np.take_along_axis(idx, part, 1)

    #if verbose:
    #    print('Compute exact, full knn')
    #ref_index = hnswlib.Index(space=metric_hnsw, dim=n_pcs)
    #ref_index.init_index(
    #    max_elements=ntot,
    #    ef_construction=ef,
    #    M=M,
    #)
    #ref_index.add_items(adatam.obsm[f'X_{space}'], np.arange(ntot))
    #ref_index.set_ef(ef)  # FIXME: Probably redundant?
    #idx, dist = ref_index.knn_query(adatam.obsm[f'X_{space}'], k=n_internal_neighbors + 1)
    ## Exclude loops/identities
    #idx = idx[:, 1:]
    #dist = dist[:, 1:]
    #dist = np.maximum(dist, 0)

    #from scipy.spatial.distance import cdist
    #idx = np.zeros((ntot, n_internal_neighbors), int)
    #dist = np.zeros((ntot, n_internal_neighbors), np.float32)
    #x = adatam.obsm[f'X_{space}']
    #i = 0
    #blk = 100
    #while i < ntot:
    #    print(f'knn vertex: {i+1} / {ntot}', end='\r')
    #    y = x[i : i + blk]
    #    d = cdist(y, x, metric=metric)

    #    # FIXME: Exclude internal neighbors for testing
    #    #d[:, n_reference:] = 1000

    #    idxi = d.argsort(axis=1)[:, 1:n_internal_neighbors + 1]
    #    disti = np.take_along_axis(d, idxi, 1)
    #    idx[i : i+blk] = idxi
    #    dist[i: i+blk] = disti
    #    i += blk

    #dist = np.maximum(dist, 0)

    ##  FIXME: no need for this I guess
    #order = dist.argsort()
    #dist = np.take_along_axis(dist, order, 1)
    #idx = np.take_along_axis(idx, order, 1)

    if verbose:
        print('Convert to sparse distance matrix')
    dist_matrix = lil_matrix((ntot, ntot), dtype=np.float32)
    for i, (js, ds) in enumerate(zip(idx, dist)):
        for j, d in zip(js, ds):
            dist_matrix[i, j] = d

    # Add incoming edges
    #for i, (js, ds) in enumerate(zip(idx_incoming, dist_incoming)):
    #    for j, d in zip(js, ds):
    #        # Hard cutoff for edges
    #        if d < 0.3:
    #            dist_matrix[i, j] = d

    if verbose:
        print('Convert lil matrix to csr')
    dist_matrix = dist_matrix.tocsr()

    if verbose:
        print('Add matrix to adata.obsp')
    adatam.obsp['distances'] = dist_matrix

    if (space == 'pca') and (features != 0):
        if verbose:
            print('Recover full gene expression matrix from .raw')
        adatam.raw.to_adata()

    return adatam
