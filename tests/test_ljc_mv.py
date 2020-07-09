# based off: https://github.com/pysal/esda/blob/master/tests/test_join_counts.py
import unittest
import numpy as np
from libpysal.weights.util import lat2W
from libpysal.common import pandas

PANDAS_EXTINCT = pandas is None

class Local_Join_Counts_MV_Tester(unittest.TestCase):
    """Unit test for Local Join Counts (univariate)"""
    def setUp(self):
        self.w = lat2W(4, 4)
        self.x = np.ones(16)
        self.x[0:8] = 0
        self.y = [0,1,0,1,1,1,1,1,0,0,1,1,0,0,1,1]
        self.z = [0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,1]


    def test_Local_Join_Counts_MV(self):
            """Test method"""
            np.random.seed(12345)
            ljc_mv = Local_Join_Count_MV(connectivity=self.w).fit([self.x, self.y, self.z])
            self.assertAlmostEqual(ljc_mv.LJC, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 2])
            
suite = unittest.TestSuite()
test_classes = [
    Local_Join_Counts_MV_Tester
]
for i in test_classes:
    a = unittest.TestLoader().loadTestsFromTestCase(i)
    suite.addTest(a)

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    runner.run(suite)