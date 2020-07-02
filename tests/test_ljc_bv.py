# based off: https://github.com/pysal/esda/blob/master/tests/test_join_counts.py
import unittest
import numpy as np
from libpysal.weights.util import lat2W
from libpysal.common import pandas

PANDAS_EXTINCT = pandas is None

class Local_Join_Counts_BV_Tester(unittest.TestCase):
    """Unit test for Local Join Counts (univariate)"""
    def setUp(self):
        self.w = lat2W(4, 4)
        self.x = np.ones(16)
        self.x[0:8] = 0
        self.z = [0,1,0,1,1,1,1,1,0,0,1,1,0,0,1,1]

    def test_Local_Join_Counts_BV(self):
            """Test method"""
            np.random.seed(12345)
            ljc_bv_case1 = Local_Join_Count_BV(connectivity=self.w).fit(self.x, self.z, case="BJC")
            ljc_bv_case2 = Local_Join_Count_BV(connectivity=self.w).fit(self.x, self.z, case="CLC")
            self.assertAlmostEqual(ljc_bv_case1.LJC, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0])
            self.assertAlmostEqual(ljc_bv_case2.LJC, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 2, 2])
            
            

