# treasuremap
![Logo](https://raw.githubusercontent.com/iosonofabio/treasuremap/master/logo.png)

Treasuremap is a variation of UMAP that embeds a network with some fixed nodes.

## Brief description
`treasuremap` is a Python package to embed a network (often a single cell omics network) derived from UMAP [1]. Unlike UMAP, you can
specify the coordinates of some nodes/cells if they are already known. They will not change, while the rest of the space warps around them.

treasuremap's superpower is using the coordinates from a known embedding (e.g. a cell atlas) and harmonizing a new data set onto that preexisting map.

Also, treasuremap was heavily developed during paternity leave after the birth of my beautiful daughter.

## Documentation
https://treasuremap.readthedocs.io

## Installation
```bash
pip install treasuremap
```

### Dependencies
`treasuremap` offers three interfaces:

1. Using `anndata`: you'll need `numpy`, and `scipy`. `scanpy` is not required but you'll probably have it anyway.
2. Using `igraph`: no additional requirements beyond the Python interface of `igraph` itself.
3. Directly setting the edges and distances of the network: in this case there are no dependencies.

`treasuremap` ships its own version of the `igraph` C core because a recent version (>=0.10.0rc2) is needed. If you use the Python `igraph` interface, that can be an older version (e.g. 0.9.x).

## Usage
Treasuremap can be used either on its own to perform only the embedding step, or it can be combined with [northstar](https://github.com/northstaratlas/northstar).

### Standalone
Standalone, `treasuremap` starts from a graph and a boolean list describing what nodes should stay fixed and computes an embedding:
```python
import igraph as ig
import treasuremap

# Generate a random similarity graph
graph = ig.Graph(n=6, edges=[
    (0, 1),
    (0, 2),
    (0, 3),
    (3, 4),
    (4, 5),
    ])
distances = [0.1, 0.3, 0.1, 0.2, 0.1]

# Initial coordinates
seed = [[-1, 1],
        [-1, 0.5],
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0]]

# Set which nodes coordinates are fixed
is_fixed = [True, True, False, False, False, False]

# Compute graph layout without touching the coordinates
# of the first two vertices
layout = treasuremap.treasuremap_igraph(
    graph, dist=distances, seed=seed, is_fixed=is_fixed,
)
```

### Together with northstar
Combining `treasuremap` with `northstar` enables a full workflow to navigate single cell data. You start from a data set that has a known cell type annotation and an existing embedding (e.g. UMAP): this is typically a cell atlas or an approximation. You then use `treasuremap` + `northstar` to analyze a new data set **onto** that template, obtaining both cell type annotation and embedding for your new data that is consistent with the atlas used.

This workflow is a little more complex, therefore `treasuremap` uses a class object called `ModelWithNorthstar`:
```python
import treasuremap

# adata with known embedding
adata_fixed = ...

# adata to co-embed
adata_new = ...

model = treasuremap.ModelWithNorthstar(
   adata_fixed, adata_new,
   seed_name='umap',
   northstar_cluster=True,
   northstar_options={'atlas_annotation_column': 'cell_type'},
)

model.build_graph()
model.cluster_graph()
model.embed_graph()

# the result is a pandas dataframe indexed by
# the cell names, first the ones from adata_fixed,
# followed by the ones from adata_new. Columns are
# one for annotation and one for each embedding dimension
result = model.result
```

## Citations
We are writing a paper describing `treasuremap`

[1] McInnes & Healy (2018), UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction, arXiv.

## License
`treasuremap` is released under the GPL 3.0 License.
Parts of the code are adapted from `igraph`, which is released under the GPL 2.0 License.
`northstar` is released under the MIT license.

