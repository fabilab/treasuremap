# vim: fdm=indent
'''
author:     Fabio Zanini
date:       17/11/22
content:    Compute a subspace of feature scores which is amenable to knn
            graph computation.
'''
import numpy as np
import pandas as pd


def compute_gene_groups(adata, cell_subtype_column, cell_type_column=None):
    '''Compute gene groups for dimensionality reduction using cell type info

    Args:
        adata: AnnData object to find gene groups for. This is typically a
          compressed cell atlas. You are not expected to call this function on
          a full dense scRNA-Seq data set.
        cell_subtype_column: A string with the name of the adata.obs column
          containing the cell subtype information.
        cell_type_column: An optional string with the name of the adata.obs
          column containing the cell type information. This level of cell
          typing is supposed to be a hierarchical parent of the
          cell_subtype_column argument. Typical cell types at this level are
          "endothelial", "epithelial", "mesenchymal", "neuron", and "immune".
          If None, only subtype information will be used, at the cost of lower
          accuracy of the dimensionality reduction.

    Returns: A dict with cell types/subtypes as keys and lists of genes as
      values. This object can be plugged into treasuremap.project_onto_subspace
      directly.
    '''
    df = adata.obs[[cell_subtype_column]].copy()
    df['idx'] = np.arange(len(df))
    if cell_type_column is not None:
        df[cell_type_column] = adata.obs[cell_type_column]

    # Fractions of cells expressing
    fractions = []

    # Fraction expressing in cell subtype
    fr_exps = {}
    for ct, idxs in df.groupby(cell_subtype_column):
        idxs = idxs['idx'].values
        fr = (adata.X[idxs] != 0).mean(axis=0)
        fr_exps[ct] = np.asarray(fr)[0]
    fr_exps = pd.DataFrame(fr_exps, index=adata.var_names)
    fractions.append(fr_exps)

    # Fraction expressing in cell type
    if cell_type_column is not None:
        fr_par = {}
        for ct, idxs in df.groupby('parent'):
            idxs = idxs['idx'].values
            fr = (adata.X[idxs] != 0).mean(axis=0)
            fr_par[ct] = np.asarray(fr)[0]
        fr_par = pd.DataFrame(fr_par, index=adata.var_names)
        fractions.append(fr_par)

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
                raise ValueError(f"No markers found for cell type {ct}")
            gene_groups[ct] = tmp.index.tolist()

    return gene_groups


def project_onto_subspace(adata, gene_groups, group_order=None, inplace=True):
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
