# based off: https://github.com/pysal/esda/blob/master/tests/test_join_counts.py
import unittest
import numpy as np
from libpysal.weights.util import lat2W
from libpysal.common import pandas

PANDAS_EXTINCT = pandas is None

class Local_Join_Counts_Tester(unittest.TestCase):
    """Unit test for Local Join Counts (univariate)"""
    def setUp(self):
        self.w = lat2W(4, 4)
        self.y = np.ones(16)
        self.y[0:8] = 0

    def test_Local_Join_Counts(self):
            """Test method"""
            np.random.seed(12345)
            ljc = Local_Join_Count(connectivity=self.w).fit(self.y)
            self.assertAlmostEqual(ljc.LJC, [0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 3, 2, 2, 3, 3, 2])