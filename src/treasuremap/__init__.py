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


def _recenter(result):
    """Recenter layout after algo"""
    if len(result) == 0:
        return
    dim = len(result[0])

    # Find center
    centers = [0 for i in range(dim)]
    for cell in result:
        for i, coord in enumerate(cell):
            centers[i] += coord
    for i in range(dim):
        centers[i] /= len(result)

    # Recenter in place
    for cell in result:
        for i, coord in enumerate(cell):
            cell[i] -= centers[i]


def _layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist=None,
        seed=None,
        is_fixed=None,
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
        negative_sampling_rate=5,
        a=-1,
        b=-1,
        distances_are_connectivities=False,
    ):

    if min_dist < 0:
        raise ValueError(f"Minimum distance must be positive, got {min_dist}")

    if epochs < 1:
        raise ValueError(f"Number of epochs must be positive, got {epochs}")

    if (sampling_prob <= 0) or (sampling_prob > 1):
        raise ValueError(f"Sampling probability must be in (0, 1], for {sampling_prob}")

    if dim not in (2, 3):
        raise ValueError(f"Number of dimensions must be 2 or 3, for {dim}")

    if nvertices == 0:
        return Layout([])
    if nvertices == 1:
        return Layout([[0, 0]])

    if len(edges) == 0:
        raise ValueError("graph has no edges, Treasuremap cannot do much...")

    if is_fixed is None:
        is_fixed = [False] * nvertices
    # If all vertices are fixed, why are you calling this?
    elif sum(is_fixed) == nvertices:
        return [list(x) for x in seed]

    if len(is_fixed) != nvertices:
        raise ValueError("is_fixed must be a boolean vector of length n. vertices")

    if dist is None:
        dist = [1] * nedges
    if len(dist) != nedges:
        raise ValueError("dist must be a vector with length n. edges")

    # Call C function
    result = _treasuremap.layout_treasuremap(
        nvertices, nedges,
        edges, dist, seed, is_fixed,
        min_dist, sampling_prob, epochs, dim,
        negative_sampling_rate,
        a, b,
        distances_are_connectivities,
    )

    return result


def treasuremap_adata(
        adata : AnnData,
        seed_name='umap',
        is_fixed=None,
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
        copy=False,
        negative_sampling_rate=5,
        seed_nonfixed='closest_fixed',
        recenter_layout=False,
    ):
    if AnnData is None:
        raise ImportError("Install the package anndata to use this function")

    nvertices = adata.shape[0]

    if 'distances' not in adata.obsp:
        raise KeyError("AnnData object must have an obsp['distances'] matrix")
    dist_matrix = adata.obsp['distances'].tocoo()
    dist_matrix.sum_duplicates()

    edges = list(zip(dist_matrix.row, dist_matrix.col))
    dist = list(dist_matrix.data)
    nedges = len(dist)
    if is_fixed is not None:
        is_fixed = [True if x else False for x in is_fixed]

    kwargs = {}
    if (seed_name is None) or (seed_name == ''):
        seed = None
    else:
        obsm_name = f'X_{seed_name}'
        if obsm_name not in adata.obsm:
            raise KeyError("AnnData object missing {obsm_name} field")
        seed = adata.obsm[f'X_{seed_name}']
        seed_dim = seed.shape[1]
        if seed_dim != dim:
            raise ValueError("{obsm_name} is {seed_dim}D, requested {dim}D embedding")
        seed = seed.tolist()

        # FIXME
        adata.obs['no_fixed_neighbor'] = False

        # If requested, seed free nodes with the coordinate of a fixed neighbor
        # We introduce a little noise to avoid exact overlapping
        if (seed_nonfixed == 'closest_fixed') and (is_fixed is not None) and any(is_fixed):
            dist_matrix = dist_matrix.tocsr()
            no_fixed_nei = []
            for k, (coords, row, fix) in enumerate(zip(seed, dist_matrix, is_fixed)):
                # If you're a fixed node, you already have fixed coords
                if fix:
                    continue
                # If you're a free node, look for neighbors that are fixed
                for i in row.indices:
                    # If a fixed neighbor is found, take its coordinates
                    if is_fixed[i]:
                        for j in range(len(coords)):
                            coords[j] = seed[i][j] + 0.01 * np.random.rand()
                        break
                else:
                    # No fixed neighbor, accumulate and seed at the end
                    no_fixed_nei.append(k)

            if len(no_fixed_nei):
                nn = len(no_fixed_nei)
                print(f'{nn} cells have no fixed neighbors')

            # FIXME
            adata.obs['no_fixed_neighbor'] = False
            adata.obs['no_fixed_neighbor'].iloc[no_fixed_nei] = True

            while no_fixed_nei:
                k = no_fixed_nei[0]
                row = dist_matrix[k]
                coords = seed[k]
                for i in row.indices:
                    # This neighbor also has no fixed neighbors, look for a
                    # better neighbor if available
                    if i in no_fixed_nei:
                        continue
                    # This neighbor has fixed neighbors, so it got assigned a
                    # coordinate already, copy it over
                    for j in range(len(coords)):
                        coords[j] = seed[i][j] + 0.01 * np.random.rand()
                    break
                else:
                    # All neighbors are also without a neighbor, use a totaly
                    # random starting position
                    for j in range(len(coords)):
                        coords[j] = 15.0 * np.random.rand()
                del no_fixed_nei[0]

        if (seed_name in adata.uns) and ('params' in adata.uns[seed_name]):
            if 'a' in adata.uns[seed_name]['params']:
                kwargs['a'] = adata.uns[seed_name]['params']['a']
            if 'b' in adata.uns[seed_name]['params']:
                kwargs['b'] = adata.uns[seed_name]['params']['a']

    result = _layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist,
        seed,
        is_fixed,
        min_dist,
        sampling_prob,
        epochs,
        dim,
        negative_sampling_rate=negative_sampling_rate,
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
        min_dist=0.01,
        sampling_prob=1.0,
        epochs=10,
        dim=2,
        negative_sampling_rate=5,
    ):

    if Graph is None:
        raise ImportError("Install the package igraph to use this function")

    nvertices = graph.vcount()
    nedges = graph.ecount()
    edges = graph.get_edgelist()

    result = _layout_treasuremap(
        nvertices,
        nedges,
        edges,
        dist,
        seed,
        is_fixed,
        min_dist,
        sampling_prob,
        epochs,
        dim,
        negative_sampling_rate=negative_sampling_rate,
    )

    # Recenter
    _recenter(result)

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
            sampling_prob=1.0,
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
            sampling_prob=sampling_prob,
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
