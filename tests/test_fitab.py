import unittest
from treasuremap._treasuremap import fit_ab


class TreasuremapTests(unittest.TestCase):
    def testFitab(self):
        # Table made with umap.umap_.fit_ab_params
        # Each row is: (min_dist, (a, b))
        expected = [
            (0.001, (1.9290733968410037, 0.7915045333956955)),
            (0.003, (1.921613194452872, 0.7935268106450729)),
            (0.01, (1.8956058664239412, 0.8006378441176886)),
            (0.03, (1.822210075237197, 0.8211992026751771)),
            (0.1, (1.5769434604035877, 0.8950608780665811)),
            (0.3, (0.9921756197688717, 1.1122533842193434)),
            (1.0, (0.11497568308423263, 1.9292371454400927)),
            ]
        for min_dist, (a_exp, b_exp) in expected:
            a, b = fit_ab(min_dist)
            self.assertEqual(a, a_exp)
            self.assertEqual(b, b_exp)
