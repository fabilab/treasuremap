import unittest
from math import hypot
import numpy as np
import scipy as sp
import anndata
from treasuremap import treasuremap_adata


class TreasuremapTests(unittest.TestCase):
    def testEmpty(self):
        adata = anndata.AnnData()
        adata.obsp['distances'] = sp.sparse.coo_matrix((0, 0), np.float32)
        adata.obsm['X_umap'] = np.zeros((0, 2), np.float32)

        self.assertRaises(
                ValueError, treasuremap_adata,
                adata=adata,
                min_dist=-0.01,
            )

        self.assertRaises(
                ValueError, treasuremap_adata,
                adata=adata,
                epochs=-1,
            )

        self.assertRaises(
                ValueError, treasuremap_adata,
                adata=adata,
                sampling_prob=-0.01,
            )

        self.assertRaises(
                ValueError, treasuremap_adata,
                adata=adata,
                sampling_prob=1.01,
            )

        self.assertRaises(
                ValueError, treasuremap_adata,
                adata=adata,
                dim=1,
            )

        # Empty anndata
        lo = treasuremap_adata(adata, copy=True)
        self.assertTrue(isinstance(lo, np.ndarray))
        self.assertEqual(lo.size, 0)

    def testSingleton(self):
        adata = anndata.AnnData(X=np.ones((1, 20), np.float32))
        adata.obsp['distances'] = sp.sparse.coo_matrix(np.ones((1, 1), np.float32))
        adata.obsm['X_umap'] = np.zeros((1, 2), np.float32)

        lo = treasuremap_adata(adata, copy=True)
        self.assertEqual(lo.tolist(), [[0, 0]])

    def testComplex(self):
        adata = anndata.AnnData(X=np.ones((12, 20), np.float32))
        adata.obsm['X_umap'] = np.zeros((12, 2), np.float32)

        # Graph with two articulation points
        edges = [
            0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3,
            3, 4, 4, 5, 5, 6,
            6, 7, 7, 8, 6, 8, 7, 9, 6, 9, 8, 9, 7, 10, 8, 10, 9, 10,
            10, 11, 9, 11, 8, 11, 7, 11,
            ]
        edges = list(zip(edges[::2], edges[1::2]))
        dist = [
            0.1, 0.09, 0.12, 0.09, 0.1, 0.1,
            0.9, 0.9, 0.9,
            0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.08, 0.05, 0.1, 0.08, 0.12, 0.09, 0.11
            ]
        dmat = sp.sparse.lil_matrix((12, 12), dtype=np.float32)
        for (i, j), dis in zip(edges, dist):
            dmat[i, j] = dis
        adata.obsp['distances'] = dmat.tocoo()

        lo = treasuremap_adata(adata, epochs=500, sampling_prob=0.3, copy=True,
                               seed_name=None)
        self.assertTrue(isinstance(lo, np.ndarray))
        print(lo)

        # One should get two clusters in this case
        x, y = list(zip(*lo))
        xmax, ymax, xmin, ymin = max(x), max(y), min(x), min(y)
        distmax = max(xmax - xmin, ymax - ymin)
        for iclu in range(0, 8, 7):
            xclu = sum(x[iclu:iclu + 4]) / 4
            yclu = sum(y[iclu:iclu + 4]) / 4
            for i in range(4):
                dx = x[iclu + i] - xclu
                dy = y[iclu + i] - yclu
                dxy = hypot(dx, dy)
                # Distance from each cluster's center should be relatively small
                self.assertLess(dxy, 0.2 * distmax)

        # Test single epoch with seed
        adata.obsm['X_umap'] = lo
        lo_adj = treasuremap_adata(adata, epochs=1, sampling_prob=1,
                                   seed_name='umap', copy=True)
        self.assertTrue(isinstance(lo_adj, np.ndarray))
