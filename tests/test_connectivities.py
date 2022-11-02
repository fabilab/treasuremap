import unittest
import numpy as np
from treasuremap._treasuremap import compute_connectivities


class TreasuremapTests(unittest.TestCase):
    def testEmpty(self):
        nvertices = 0
        edges = []
        distances = []
        connectivities = compute_connectivities(
                nvertices, len(edges), edges, distances)

        self.assertEqual(type(connectivities), list)
        self.assertEqual(len(connectivities), 0)

    def testNullDistance(self):
        nvertices = 2
        edges = [(0, 1)]
        connectivities = compute_connectivities(
                nvertices, len(edges), edges, None)

        self.assertEqual(type(connectivities), list)
        self.assertEqual(len(connectivities), 1)
        self.assertEqual(connectivities[0], 1.0)

    def testSingleton(self):
        nvertices = 2
        edges = [(0, 1)]
        distances = [14.0]
        connectivities = compute_connectivities(
                nvertices, len(edges), edges, distances)

        self.assertEqual(type(connectivities), list)
        self.assertEqual(len(connectivities), 1)
        self.assertEqual(connectivities[0], 1.0)
