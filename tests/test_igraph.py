import unittest
from math import hypot
import numpy as np
from igraph import Graph, Layout, InternalError
from treasuremap import treasuremap_igraph


class TreasuremapTests(unittest.TestCase):
    def testEmpty(self):
        g = Graph()

        self.assertRaises(
                ValueError, treasuremap_igraph,
                graph=g,
                min_dist=-0.01,
            )

        self.assertRaises(
                ValueError, treasuremap_igraph,
                graph=g,
                epochs=-1,
            )

        self.assertRaises(
                ValueError, treasuremap_igraph,
                graph=g,
                negative_sampling_rate=-1,
            )

        self.assertRaises(
                ValueError, treasuremap_igraph,
                graph=g,
                dim=1,
            )

        # Empty graph
        lo = treasuremap_igraph(g)
        self.assertTrue(isinstance(lo, Layout))
        self.assertEqual(lo.coords, [])

    def testSingleton(self):
        # Singleton graph
        g = Graph(n=1)
        lo = treasuremap_igraph(g, recenter_layout=True)
        self.assertEqual(lo.coords, [[0, 0]])

    def testSmall(self):
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
            0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.08, 0.01, 0.1, 0.08, 0.12, 0.09, 0.11
            ]
        g = Graph(edges)
        lo = treasuremap_igraph(g, dist=dist, epochs=500)
        self.assertTrue(isinstance(lo, Layout))

        # One should get two clusters in this case
        x, y = list(zip(*lo.coords))
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
                self.assertLess(dxy, 0.4 * distmax)

        # Test single epoch with seed
        lo_adj = treasuremap_igraph(g, dist=dist, epochs=1, seed=lo)
        self.assertTrue(isinstance(lo_adj, Layout))

        # Same but inputting the coordinates
        lo_adj = treasuremap_igraph(g, dist=dist, epochs=1, seed=lo.coords)
        self.assertTrue(isinstance(lo_adj, Layout))

    def testLarge(self):
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

        g = Graph(edges)
        lo = treasuremap_igraph(
            g,
            dist=distances,
            epochs=200,
        )
        self.assertTrue(isinstance(lo, Layout))

        n_vertices = 0
        for cloud_size in cloud_sizes:
             err = np.std(lo[n_vertices: n_vertices + cloud_size], axis=0)
             self.assertTrue((err < 1.0).all())
             n_vertices += cloud_size
