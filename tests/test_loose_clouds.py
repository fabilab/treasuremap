# vim: fdm=indent
'''
author:     Fabio Zanini
date:       31/08/22
content:    Test algo on fly atlas data.
'''
import sys
import time
import numpy as np
import pandas as pd
import igraph as ig
import treasuremap

import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':

    # Loosely connected clouds of highly connected vertices
    cloud_sizes = [50, 100, 30, 20]
    n_vertices = 0
    edges = []
    distances = []
    for cloud_size in cloud_sizes:
        for v1 in range(cloud_size):
            for v2 in range(1, 11):
                edges.append(
                    (v1 + n_vertices, v2 + n_vertices),
                )
                distances.append(2)
        n_vertices += cloud_size

    g = ig.Graph(edges)
    lo = treasuremap.treasuremap_igraph(
        g,
        dist=distances,
        epochs=50,
    )

    fig, ax = plt.subplots()
    ig.plot(g, target=ax, layout=lo)
    fig.tight_layout()

    coords = np.array(lo.coords)

    clouds = pd.Series(np.zeros(n_vertices, np.int32))
    clouds[100:250] = 1
    clouds[250:280] = 2
    clouds[280:] = 3
    clouds = clouds.to_frame(name='cloud')
    clouds['c'] = 1

    print('Plot embeddings')
    cmap = dict(zip(np.arange(4), sns.color_palette('husl', n_colors=4)))
    fig, ax = plt.subplots(figsize=(4, 4))
    for cell_type, group in clouds.groupby('cloud'):
        cellnames = group.index
        x, y = coords[cellnames].T
        ax.scatter(x, y, color=cmap[cell_type], alpha=0.5, label=str(cell_type))

    ax.set_title('Treasuremap (free)')
    ax.legend()
    fig.tight_layout()
    plt.ion(); plt.show()
