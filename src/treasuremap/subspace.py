# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/11/22
content:    Compute a subspace of feature scores which is amenable to knn
            graph computation.
'''
import numpy as np
import pandas as pd


def compute_gene_groups(adata, cell_type_columns):
    '''Compute gene groups for dimensionality reduction using cell type info

    Args:
        adata: AnnData object to find gene groups for. This is typically a
          compressed cell atlas. You are not expected to call this function on
          a full dense scRNA-Seq data set.
        cell_type_columns: A column name or list of column names that are in
          adata.obs, to be used to compute gene groups. These are typically
          the columns containing the cell type and subtype annotations, e.g.
          cell_type_columns=['celltype', 'cellsubtype'].

    Returns: A dict with cell types/subtypes as keys and lists of genes as
      values. This object can be plugged into treasuremap.project_onto_subspace
      directly.
    '''
    if isinstance(cell_type_columns, str):
        cell_type_columns = [cell_type_columns]
    df = adata.obs[cell_type_columns].copy()
    df['idx'] = np.arange(len(df))

    # Fractions of cells expressing
    fractions = []

    # Fraction expressing in each column (e.g. within a cell type)
    for column in cell_type_columns:
        fr_exps = {}
        for ct, idxs in df.groupby(column):
            idxs = idxs['idx'].values
            fr = (adata.X[idxs] != 0).mean(axis=0)
            fr_exps[ct] = np.asarray(fr)[0]
        fr_exps = pd.DataFrame(fr_exps, index=adata.var_names)
        fractions.append(fr_exps)

    # Find markers for each cell subtype/type, and add them up until a certain
    # KS statistic. This might change in the future, hence kssum is not public
    kssum = 5
    gene_groups = {}
    for fr in fractions:
        for ct in fr.columns:
            other = [col for col in fr.columns if col != ct]
            diff = -fr[other].copy()
            for col in other:
                diff[col] += fr[ct].values

            # Get the most similar type, then the genes with the largest distance
            # to the closest type
            tmp = diff.min(axis=1).nlargest(200).cumsum()
            tmp = tmp[:(tmp >= kssum).values.nonzero()[0][0] + 1]
            if len(tmp) == 0:
                raise ValueError(f"No markers found for {ct}")
            gene_groups[ct] = tmp.index.tolist()

    return gene_groups


def project_onto_subspace(
        adata,
        gene_groups,
        group_order=None,
        inplace=True,
        normalize_by_row=False,
        ):
    '''Compute subspace from gene groups

    Args:
        adata: AnnData object to compute a subspace for
        gene_groups: a dict with cell types/subtypes as keys and gene lists
          as values. Each dimension of the subspace will correspond to one
          gene group.
        group_order: a list of cell types/subtypes that specifies the order
          of the gene groups in the subspace matrix. If None, the default
          Python order for dictionaries is used, and the actual order is
          returned/stored.
        inplace: Whether the adata object is modified in place. If True,
          adata.obsm['X_subspace'] will contain the matrix and
          adata.unsw['subspace'] the cell type names and gene groups. If
          False, a dictionary with those things is returned.
    '''


    x = adata.X.tocsc()
    is_log = x.max() < 50
    if not is_log:
        x.data = np.log(x.data + 1)

    #x_avg = x.mean(axis=0)
    #x_std = np.sqrt(np.power((x - x_avg), 2).sum(axis=0))

    n, m = adata.shape
    ndim = len(gene_groups)
    mat = np.zeros((n, ndim), np.float32)

    genes = pd.Series(np.arange(m), index=adata.var_names)
    if group_order is None:
        group_order = list(gene_groups.keys())
    for i, gn in enumerate(group_order):
        mki = gene_groups[gn]
        idxi = genes.loc[mki].values
        xi = x[:, idxi].copy()

        # Cap from above
        #xi.data = np.minimum(xi.data, 2)

        # Convert to standard scale
        xi /= (xi.max(axis=0).toarray()[0] + 1e-3)

        # Convert to zscale
        #xi = (xi - x_avg[:, idxi]) / (x_std[:, idxi] + 1e-3)

        scorei = np.asarray(xi.sum(axis=1))[:, 0]
        mat[:, i] = scorei

    if normalize_by_row:
        mat /= (mat.sum(axis=1) + 1e-3)

    if inplace:
        adata.obsm['X_subspace'] = mat
        adata.uns['subspace'] = {
            'group_order': group_order,
            'gene_groups': gene_groups,
        }
    else:
        return {
            'group_order': group_order,
            'gene_groups': gene_groups,
            'X': mat,
        }
