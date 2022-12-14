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
from treasuremap.compress_atlas import subsample_atlas
from treasuremap.knn import build_knn
from treasuremap.weights import compute_weights
from treasuremap.subspace import compute_gene_groups, project_onto_subspace


def treasuremap_adata(
        adata : AnnData,
        seed_name='umap',
        is_fixed=None,
        min_dist=0,
        epochs=100,
        dim=2,
        copy=False,
        negative_sampling_rate=4,
        seed_nonfixed='closest_fixed',
        recenter_layout=False,
        use_weights=False,
    ):
    """Compute a Treasuremap embedding starting from an AnnData object

    Args:
    ...
    """

    if AnnData is None:
        raise ImportError("Install the package anndata to use this function")

    # NOTE: AnnData does not store the graph per se, only the distances and
    # usually the connectivities as weighted, sparse matrices. So there is no
    # obvious way to use treasuremap's ability to work on unweighted graphs
    # when using the AnnData interface.
    if not use_weights:
        if 'distances' not in adata.obsp:
            raise KeyError("AnnData object must have an obsp['distances'] matrix")
        dist_matrix = adata.obsp['distances'].tocoo()
        # NOTE: dist_matrix is basically a directed graph, i.e. it is not symmetric
        # Later, when it is converted into weights, it gets symmetrized. Nonetheless
        # it is conceivable that the technical sparse matrix has duplicates that
        # would result in parallel edges. We should get rid of them.
        dist_matrix.sum_duplicates()
    else:
        if 'connectivities' not in adata.obsp:
            raise KeyError("AnnData object must have an obsp['connectivities'] matrix")
        dist_matrix = adata.obsp['connectivities'].tocoo()

    # Get edges and distances, ignoring loops. The edges are still basically
    # directed, so no symmetry yet.
    edges = []
    dist = []
    for i, j, d in zip(dist_matrix.row, dist_matrix.col, dist_matrix.data):
        if i == j:
            continue
        edges.append((i, j))
        dist.append(d)
    nedges = len(dist)
    nvertices = adata.shape[0]

    if is_fixed is not None:
        is_fixed = [True if x else False for x in is_fixed]

    kwargs = {}
    if (seed_name is None) or (seed_name == ''):
        seed = np.random.rand(nvertices, dim).tolist()
    else:
        obsm_name = f'X_{seed_name}'
        if obsm_name not in adata.obsm:
            raise KeyError(f"AnnData object missing {obsm_name} field")
        seed = adata.obsm[f'X_{seed_name}']
        seed_dim = seed.shape[1]
        if seed_dim != dim:
            raise ValueError("{obsm_name} is {seed_dim}D, requested {dim}D embedding")
        seed = seed.tolist()

        # If requested, seed free nodes with the coordinate of a fixed neighbor
        # If the vertex has no fixed neighbor, just put it at random
        if (seed_nonfixed == 'closest_fixed') and (is_fixed is not None) and any(is_fixed):
            dist_matrix = dist_matrix.tocsr()
            for k, (coords, row, fix) in enumerate(zip(seed, dist_matrix, is_fixed)):
                # If you're a fixed node, you already have fixed coords
                if fix:
                    continue
                # If you're a free node, look for neighbors that are fixed
                for i in row.indices:
                    # If a fixed neighbor is found, take its coordinates
                    if is_fixed[i]:
                        seed[k] = [seed[i][j] for j in range(dim)]
                        break
                else:
                    seed[k] = np.random.rand(dim).tolist()

        if (seed_name in adata.uns) and ('params' in adata.uns[seed_name]):
            if 'a' in adata.uns[seed_name]['params']:
                kwargs['a'] = adata.uns[seed_name]['params']['a']
            if 'b' in adata.uns[seed_name]['params']:
                kwargs['b'] = adata.uns[seed_name]['params']['b']

    result = _treasuremap.layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist,
        seed,
        is_fixed,
        min_dist,
        epochs,
        dim,
        negative_sampling_rate=negative_sampling_rate,
        distances_are_weights=use_weights,
        **kwargs,
    )
    result = np.asarray(result).astype(np.float32)

    if recenter_layout and len(result):
        result -= result.mean(axis=0)

    if copy:
        return result

    adata.obsm['X_treasuremap'] = result


def treasuremap_igraph(
        graph : Graph,
        dist=None,
        seed=None,
        is_fixed=None,
        min_dist=0,
        epochs=100,
        dim=2,
        negative_sampling_rate=5,
        recenter_layout=False,
    ):
    """Compute a Treasuremap embedding starting from an igraph Graph."""

    if Graph is None:
        raise ImportError("Install the package igraph to use this function")

    # Make sure there are no loops or parallel edges. This requires copies so
    # we only do it if necessary
    if any(graph.is_loop()) or graph.has_multiple():
        graph = graph.copy()
        graph.es['treasuremap_dist'] = dist
        graph = graph.simplify(loops=True, combine_edges='min')
        dist = graph.es['treasuremap_dist']

    nvertices = graph.vcount()
    nedges = graph.ecount()
    edges = graph.get_edgelist()
    seed = np.random.rand(nvertices, dim).tolist()

    result = _treasuremap.layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist,
        seed,
        is_fixed,
        min_dist,
        epochs,
        dim,
        negative_sampling_rate=negative_sampling_rate,
    )
    result = np.asarray(result).astype(np.float32)

    if recenter_layout and len(result):
        result -= result.mean(axis=0)

    result = Layout(result)
    return result



class ModelWithNorthstar:
    model = None
    result = None

    def __init__(
        self,
        adata_fixed, adata_new,
        northstar_options=None,
        ):
        self.adata_fixed = adata_fixed
        self.adata_new = adata_new
        self.northstar_options = northstar_options

    def build_graph(self):
        import northstar

        self.model = northstar.Subsample(
            self.adata_fixed,
            **self.northstar_options,
        )

        # Run first part of northstar
        self.model.new_data = self.adata_new
        self.model._check_init_arguments()
        self.model.compute_feature_intersection()
        self.model._check_feature_intersection()
        self.model.prepare_feature_selection()
        self.model.select_features()
        self.model._check_feature_selection()
        self.model.merge_atlas_newdata()
        self.model.compute_pca()
        self.model.compute_similarity_graph()

    def cluster_graph(self):
        if self.model is None:
            raise RuntimeError("You must call build_graph first")

        self.model.cluster_graph()

        if self.result is None:
            self.result = self.model.adata_merge.obs[[self.model.atlas_annotation_column]]

    def embed_graph(
            self,
            seed_name='umap',
            min_dist=0.01,
            epochs=10,
        ):
        import pandas as pd

        if self.model is None:
            raise RuntimeError("You must call build_graph first")

        adata_fixed = self.adata_fixed
        adata_new = self.adata_new
        seed_name

        # Prepare seed
        obsm_name = f'X_{seed_name}'
        if obsm_name not in adata_fixed.obsm:
            raise KeyError(f'{obsm_name} not in adata_fixed.obsm')
        n_fixed = adata_fixed.shape[0]
        dim = adata_fixed.obsm[obsm_name].shape[1]
        seed = np.zeros((n_fixed+ adata_new.shape[0], dim), np.float32)
        seed[:n_fixed] = adata_fixed.obsm[obsm_name]
        seed = seed.tolist()

        # Prepare is_fixed
        is_fixed = np.zeros(seed.shape[0], bool)
        is_fixed[:n_fixed] = True
        is_fixed = is_fixed.tolist()

        # Get graph from northstar
        graph = self.model.graph
        if 'distance' in graph.edge_attributes:
            dist = graph.es['distance']
        else:
            dist = None

        # Embed graph
        layout = treasuremap_igraph(
            graph,
            dist=dist,
            seed=seed,
            is_fixed=is_fixed,
            min_dist=min_dist,
            epochs=epochs,
            dim=dim,
        )

        index = adata_fixed.obs_names.tolist() + adata_new.obs_names.tolist()
        columns = ['treasuremap_'+str(i+1) for i in range(dim)]
        result = pd.DataFrame(
            layout.coords, columns=columns, index=index,
        )

        if self.result is None:
            self.result = result
        else:
            for col in columns:
                self.result[col] = result[col]


__all__ = (
    'treasuremap_igraph',
    'treasuremap_adata',
    'ModelWithNorthstar',
    'subsample_atlas',
)
