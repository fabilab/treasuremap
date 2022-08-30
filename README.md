# treasuremap
Treasuremap is a variation of UMAP that embeds a network with some fixed nodes.

## Brief description
`treasuremap` is a Python package to embed a network (often a single cell omics network) similar to UMAP [1]. Unlike UMAP, you can
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
The typical u
```python
import northstar

# Choose an atlas
atlas_name = 'Darmanis_2015'

# Get a gene expression matrix of the new dataset (here a
# random matrix for simplicity)
N = 200
L = 50
new_dataset = pd.DataFrame(
    data=np.random.rand(L, N).astype(np.float32),
    index=<gene_list>,
    columns=['cell_'+str(i+1) for i in range(N)],
    )

# Initialize northstar classes
model = northstar.Averages(
        atlas='Darmanis_2015',
        n_neighbors=5,
        n_pcs=10,
        )

# Run the classifier
model.fit(new_dataset)

# Get the cluster memberships for the new cells
membership = model.membership
```

## Citation
If you use this software please cite the following paper:

Fabio Zanini\*, Bojk A. Berghuis\*, Robert C. Jones, Benedetta Nicolis di Robilant, Rachel Yuan Nong, Jeffrey Norton, Michael F. Clarke, Stephen R. Quake. **Northstar enables automatic classification of known and novel cell types from tumor samples.** Scientific Reports 10, Article number: 15251 (2020), DOI: https://doi.org/10.1038/s41598-020-71805-1

## License
`northstar` is released under the MIT license.

NOTE: The module leidenalg to perform graph-based clstering is released
under the GLP3 license. You agree with those licensing terms if you use
leidenalg within northstar.
