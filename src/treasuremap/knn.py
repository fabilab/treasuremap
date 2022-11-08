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
    features=1000,
    n_pcs=50,
    n_internal_neighbors=14,
    n_external_neighbors=5,
    verbose=0,
    metric='euclidean',
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

    if verbose:
        print('Select features')
    if isinstance(features, int):
        sc.pp.highly_variable_genes(adatam, n_top_genes=features)
    else:
        adatam.var['highly_variable'] = False
        adatam.var.loc[features, 'highly_variable'] = True

    adatam.raw = adatam
    adatam = adatam[:, adatam.var.highly_variable]

    if verbose:
        print('Scale gene expression to 10')
    sc.pp.scale(adata, max_value=10)

    if verbose:
        print('Principle Component Analysis')
    sc.tl.pca(adatam, n_comps=n_pcs)

    if verbose:
        print('Build knn')
    # Build knn in two steps: internal neighbors and external neighbors
    ef, M = 200, 48  # FIXME: check this
    metric_hnsw = 'cosine' if metric == 'correlation' else 'l2'

    #if verbose:
    #    print('Internal neighbors')
    #ref_index = hnswlib.Index(space=metric_hnsw, dim=n_pcs)
    #ref_index.init_index(
    #    max_elements=n,
    #    ef_construction=ef,
    #    M=M,
    #)
    #ref_index.add_items(adatam.obsm['X_pca'][n_reference:], np.arange(n) + n_reference)
    #ref_index.set_ef(ef)  # FIXME: Probably redundant?
    #idx_internal, dist_internal = ref_index.knn_query(
    #        adatam.obsm['X_pca'][n_reference:], k=n_internal_neighbors + 1)
    ## Exclude loops/identities
    #idx_internal = idx_internal[:, 1:]
    #dist_internal = dist_internal[:, 1:]

    ## Sometimes the feature selection results in zero expression, 
    #dist_internal = np.maximum(dist_internal, 0)

    ## Incoming edges from reference
    #idx_incoming, dist_incoming = ref_index.knn_query(
    #        adatam.obsm['X_pca'][:n_reference], k=n_external_neighbors)
    #dist_incoming = np.maximum(dist_incoming, 0)

    #if verbose:
    #    print('External neighbors')
    #ref_index = hnswlib.Index(space=metric_hnsw, dim=n_pcs)
    #ref_index.init_index(
    #    max_elements=n_reference,
    #    ef_construction=ef,
    #    M=M,
    #)
    #ref_index.add_items(adata.obsm['X_pca'][:n_reference], np.arange(n_reference))
    #ref_index.set_ef(ef)  # FIXME: Probably redundant?
    #idx_external, dist_external = ref_index.knn_query(
    #        adatam.obsm['X_pca'][n_reference:], k=n_external_neighbors)
    #dist_external = np.maximum(dist_external, 0)

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

    if verbose:
        print('Compute full knn')
    #ref_index = hnswlib.Index(space=metric_hnsw, dim=n_pcs)
    #ref_index.init_index(
    #    max_elements=ntot,
    #    ef_construction=ef,
    #    M=M,
    #)
    #ref_index.add_items(adatam.obsm['X_pca'], np.arange(ntot))
    #ref_index.set_ef(ef)  # FIXME: Probably redundant?
    #idx, dist = ref_index.knn_query(adatam.obsm['X_pca'], k=n_internal_neighbors + 1)
    ## Exclude loops/identities
    #idx = idx[:, 1:]
    #dist = dist[:, 1:]
    #dist = np.maximum(dist, 0)

    from scipy.spatial.distance import cdist
    idx = np.zeros((ntot, n_internal_neighbors), int)
    dist = np.zeros((ntot, n_internal_neighbors), np.float32)
    x = adatam.obsm['X_pca']
    i = 0
    blk = 200
    while i < ntot:
        print(f'knn vertex: {i+1} / {ntot}', end='\r')
        y = x[i : i + blk]
        d = cdist(y, x, metric=metric)

        idxi = d.argsort(axis=1)[:, 1:n_internal_neighbors + 1]
        disti = np.take_along_axis(d, idxi, 1)
        idx[i : i+blk] = idxi
        dist[i: i+blk] = disti
        i += blk

    dist = np.maximum(dist, 0)

    ##  FIXME: no need for this I guess
    #order = dist.argsort()
    #dist = np.take_along_axis(dist, order, 1)
    #idx = np.take_along_axis(idx, order, 1)

    if verbose:
        print('Convert to sparse distance matrix')
    dist_matrix = lil_matrix((ntot, ntot), dtype=np.float32)
    #for i, (js, ds) in enumerate(zip(idx, dist), n_reference):
    for i, (js, ds) in enumerate(zip(idx, dist)):
        for j, d in zip(js, ds):
            dist_matrix[i, j] = d

    # Add incomig edges
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

    if verbose:
        print('Recover full gene expression matrix from .raw')
    adatam.raw.to_adata()

    return adatam
